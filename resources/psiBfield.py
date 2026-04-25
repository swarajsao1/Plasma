import os 

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

file_path  = os.path.join(BASE_DIR, "g_file.in")

from freeqdsk import geqdsk

with open(file_path, 'r') as f:
    g = geqdsk.read(f)


# R and Z Coordinates
r_coords = g.r_grid[:,0]
z_coords = g.z_grid[0,:]

psi_norm = (g.psi - g.simagx)/(g.psi.max() - g.psi.min())


# Radial and Vertical Field Calculation
import numpy as np
dpsi_dR = np.gradient(psi_norm, r_coords, axis = 0)
dpsi_dZ = np.gradient(psi_norm, z_coords, axis =1)

BR = - (1/r_coords) * dpsi_dZ
BZ = (1/r_coords) * dpsi_dR


# Poloidal Field Calculation
from scipy.interpolate import interp1d
psi_grid = np.linspace(g.psi.min(), g.psi.max(), len(g.fpol))
fpol_interp = interp1d(psi_grid, g.fpol)

Bphi = fpol_interp(g.psi)/r_coords

print("g_file keywords extracted, BR Bphi BZ calculated.")