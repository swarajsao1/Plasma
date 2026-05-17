import os 
import sys
import numpy as np
from freeqdsk import geqdsk
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator


def extract_g(file_name="g_file.in"):

    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(BASE_DIR, file_name)

    with open(file_path, 'r') as f:
        g = geqdsk.read(f)

    print("g_file keywords extracted")
    return g


def compute_B(g):

    # Coordinates
    r_coords = g.r_grid[:, 0]
    z_coords = g.z_grid[0, :]

    # Normalized psi (for plotting only)
    psi_norm = (g.psi - g.simagx) / (g.sibdry - g.simagx)

    # Gradients
    dpsi_dR = np.gradient(g.psi, r_coords, axis=0, edge_order=2)
    dpsi_dZ = np.gradient(g.psi, z_coords, axis=1, edge_order=2)

    # Fields
    BR = -(1 / g.r_grid) * dpsi_dZ
    BZ = (1 / g.r_grid) * dpsi_dR

    # Toroidal field
    psi_grid = np.linspace(g.psi.min(), g.psi.max(), len(g.fpol))
    fpol_interp = interp1d(psi_grid, g.fpol)

    Bphi = fpol_interp(g.psi) / g.r_grid

    

    print("Psi Normalized and Magnetic field BR Bphi BZ computed.")

    return psi_norm, BR, Bphi, BZ



def field_tracing(g, psi_norm, BR, BZ):

    # =========================
    # USER INPUTS
    # =========================

    r_start = float(input("Enter r_start : ") or g.rmaxis + 0.1)
    r_end   = float(input("Enter r_end   : ") or g.rbdry.max() - 0.1)
    N       = int(input("Enter number of radial points N : ") or 10)

    initial_h = float(
        input("Enter RK4 step size initial_h [default=-0.001] : ") or -0.001)

    tolerance = float(
        input("Enter psi correction tolerance [default=1e-2] : ") or 1e-2)

    max_step = int(
        input("Enter maximum RK4 steps [default=60000] : ") or 60000)


    # ===================
    # Interpolators
    # ===================

    BR_interp = RegularGridInterpolator((g.r_grid[:,0], g.z_grid[0,:]), BR)
    BZ_interp = RegularGridInterpolator((g.r_grid[:,0], g.z_grid[0,:]), BZ)
    psi_interp = RegularGridInterpolator((g.r_grid[:,0], g.z_grid[0,:]), psi_norm)


    # =========================
    # Field Line Equation
    #==========================

    def F(r_current, z_current):
        if not (g.rbdry.min() <= r_current <= g.rbdry.max() and
            g.zbdry.min() <= z_current <= g.zbdry.max()):
            return np.nan, np.nan
        
        BR_local = BR_interp((r_current, z_current))
        BZ_local = BZ_interp((r_current, z_current))
        Bp_local = np.sqrt(BR_local**2 + BZ_local**2)
        if Bp_local == 0.0:
            br = 0.0
            bz = 0.0
        else:
            br = BR_local/Bp_local
            bz = BZ_local/Bp_local
        return br, bz


    #=========================
    # RK4 Implementation
    #=========================

    def RK4(F, r_current, z_current, h):
        k1r, k1z = F(r_current, z_current)
        if np.isnan(k1r) or np.isnan(k1z): return np.nan, np.nan

        k2r, k2z = F(r_current + 0.5*h*k1r, z_current + 0.5*h*k1z)
        if np.isnan(k2r) or np.isnan(k2z): return np.nan, np.nan

        k3r, k3z = F(r_current + 0.5*h*k2r, z_current + 0.5*h*k2z)
        if np.isnan(k3r) or np.isnan(k3z): return np.nan, np.nan

        k4r, k4z = F(r_current + h*k3r, z_current + h*k3z)
        if np.isnan(k4r) or np.isnan(k4z): return np.nan, np.nan

        r_new = r_current + (h/6)*(k1r + 2*k2r + 2*k3r + k4r)
        z_new = z_current + (h/6)*(k1z + 2*k2z + 2*k3z + k4z)
        return r_new, z_new


    # =========================
    # initial cinditions
    # =========================

    r = np.linspace(r_start, r_end, N)
    z_start = g.zmaxis
   

    r_trajectory = []
    z_trajectory = []

    # =========================
    # Field Line Tracing
    # =========================

    for i, r_initial in enumerate(r):
        r_current = r_initial
        z_current = z_start

        psi_current = psi_interp((r_current, z_current))

        r_trajectory.append([r_current])
        z_trajectory.append([z_current])

        step_count = 0
        diff_angle = 0.0


        while (step_count < max_step):

            step_count += 1
            h_current_step = initial_h

            angle_initial = np.arctan2(z_current - g.zmaxis, r_current - g.rmaxis) + np.pi
            refined_step_found = False
            for _ in range(15):
                r_trial, z_trial = RK4(F, r_current, z_current, h_current_step)

                psi_trial = psi_interp((r_trial, z_trial))

                if abs(psi_trial - psi_current) < tolerance:
                    refined_step_found = True
                    r_new, z_new = r_trial, z_trial
                    break
                else:
                    h_current_step /= 2.0
                    if h_current_step < 1e-7:
                        refined_step_found = True
                        r_new, z_new = r_trial, z_trial
                        break

            angle_trial = np.arctan2(z_trial - g.zmaxis, r_trial - g.rmaxis) + np.pi

            diff_angle += abs(angle_initial - angle_trial)

            if diff_angle >= 4*np.pi:
                break

            r_current, z_current = r_new, z_new

            r_trajectory[i].append(r_new)
            z_trajectory[i].append(z_new)


    print("Poloidal field line trajectories computed.")
    return r_trajectory, z_trajectory




def grid_points(r_trajectory, z_trajectory):
    
    cumulative_distances = []
    for i in range(len(r_trajectory)):
        dr = np.diff(r_trajectory[i])
        dz = np.diff(z_trajectory[i])
        distances = np.sqrt(dr**2 + dz**2)
        cumulative_distance = np.cumsum(distances)
        cumulative_distance = np.insert(cumulative_distance, 0, 0)
        cumulative_distances.append(cumulative_distance)

    print("Cumulative distances along trajectories computed.")
    
    N0 = int(input("Enter number of grid points along first trajectory N_grid [default=10] : ") or 10)
    
    tol = float(input("Enter grid spacing tolerance tol [default=0.5] : ") or 0.5)
    
    total_length_0 = cumulative_distances[0][-1]

    grid_distance_0 = np.linspace(0, total_length_0, N0)
    
    ds0 = grid_distance_0[1] - grid_distance_0[0]
    
    r_grid_0 = np.interp(grid_distance_0, cumulative_distances[0], r_trajectory[0])
    z_grid_0 = np.interp(grid_distance_0, cumulative_distances[0], z_trajectory[0])

    
    r_grid_points = [r_grid_0]
    z_grid_points = [z_grid_0]


    for i in range(1, len(r_trajectory)):

        cumulative_distance = cumulative_distances[i]

        total_length = cumulative_distance[-1]

        N = int(total_length / ds0) + 1

        grid_distances = np.linspace(0, total_length, N)

        ds = grid_distances[1] - grid_distances[0]

        while abs(ds - ds0) / ds0 > tol:

            N += 1

            grid_distances = np.linspace(0, total_length, N)

            ds = grid_distances[1] - grid_distances[0]


        r_grid = np.interp(grid_distances, cumulative_distance, r_trajectory[i])
        z_grid = np.interp(grid_distances, cumulative_distance, z_trajectory[i])
        
        r_grid_points.append(r_grid)
        z_grid_points.append(z_grid)

    print("Grid points along trajectories computed.")
    return r_grid_points, z_grid_points


    