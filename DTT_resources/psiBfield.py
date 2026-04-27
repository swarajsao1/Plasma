import os 
import sys
import numpy as np
from freeqdsk import geqdsk
from scipy.interpolate import interp1d


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
    psi_norm = (g.psi - g.simagx) / (g.psi.max() - g.psi.min())

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