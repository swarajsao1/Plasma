#!/usr/bin/env python
# coding: utf-8

# In[1]:


# from freeqdsk import geqdsk

# with open("g_file.in", 'r') as f:
#     g = geqdsk.read(f)


# In[2]:


from resources import psiBfield
import importlib
importlib.reload(psiBfield)

g = psiBfield.extract_g()


# In[3]:


g.z_grid[0, :]


# In[4]:


import matplotlib.pyplot as plt

plt.figure(figsize=(5,5))
plt.plot(g.rmaxis, g.zmaxis, 'r', marker = 'o', markersize = 5, label = "Magnetic Axis")
plt.plot(g.rcentr, g.zmid, 'b',marker = '^', markersize = 5, label = "Mid Plane")
plt.plot(g.rbdry, g.zbdry, 'y--', label = "Separatrix")
plt.plot(g.rlim, g.zlim, 'm', label = "Limiter")
plt.title("Tokamak Poloidal Section")
plt.xlim(g.r_grid[:,0][0], g.r_grid[:,0][-1])
plt.ylim(g.z_grid[0,:][0], g.z_grid[0,:][-1])
plt.xlabel("R [m]")
plt.ylabel("Z [m]")
plt.axis('equal')
plt.legend()
plt.show()


# In[5]:


g.rdim, g.nr


# In[6]:


# import numpy as np
# r_coords = np.linspace(g.rleft, g.rleft + g.rdim, g.nr)
# z_coords = np.linsapce(g.zmid - g.zdim/2 + g.zmid + g.zdim/2, g.nz)


# In[7]:


r_coords = g.r_grid[:,0]
z_coords = g.z_grid[0,:]


# In[8]:


plt.figure(figsize = (5,5))

plt.contour(g.r_grid, g.z_grid, g.psi, levels = 30, cmap = 'plasma')
plt.colorbar(label = r'$\psi$')
plt.xlabel("R [m]")
plt.ylabel("Z [m]")
plt.title(r"Poloidal Flux $\psi$")
plt.axis('equal')
plt.legend()
plt.show()


# In[9]:


psi_norm = (g.psi - g.simagx)/(g.psi.max() - g.psi.min())


# In[10]:


plt.figure(figsize = (5,5))

plt.contour(g.r_grid, g.z_grid, psi_norm, levels = 30, cmap = 'plasma')
plt.colorbar(label = r'$\psi$')
plt.xlabel("R [m]")
plt.ylabel("Z [m]")
plt.title(r"Normalized Poloidal Flux $\psi$")
plt.axis('equal')
plt.legend()
plt.show()


# In[11]:


import numpy as np
dpsi_dR = np.gradient(psi_norm, r_coords, axis = 0)
dpsi_dZ = np.gradient(psi_norm, z_coords, axis =1)


# In[12]:


BR = - (1/r_coords) * dpsi_dZ
BZ = (1/r_coords) * dpsi_dR


# In[13]:


plt.figure(figsize=(8,4))

plt.subplot(1,2,1)
plt.contourf(g.r_grid, g.z_grid, BR, cmap = 'plasma')
plt.colorbar(label = r"$B_R$")
plt.xlabel('R [m]')
plt.ylabel('Z [m]')
plt.axis('equal')
plt.title(r'Contour Plot of $B_R$')

plt.subplot(1,2,2)
plt.contourf(g.r_grid, g.z_grid, BZ, cmap = 'plasma')
plt.colorbar(label = r"$B_Z$")
plt.xlabel('R [m]')
plt.ylabel('Z [m]')
plt.axis('equal')
plt.title(r'Contour Plot of $B_Z$')

plt.legend()
plt.tight_layout()
plt.show()


# In[14]:


from scipy.interpolate import interp1d

psi_grid = np.linspace(g.psi.min(), g.psi.max(), len(g.fpol))
fpol_interp = interp1d(psi_grid, g.fpol)

B_phi = fpol_interp(g.psi)/r_coords


# In[15]:


plt.figure(figsize = (12,4))

plt.subplot(1,3,1)
plt.contourf(g.r_grid, g.z_grid, BR, cmap = 'plasma')
plt.colorbar(label = r'$B_R$')
plt.xlabel('R [m]')
plt.ylabel('Z [m]')
plt.axis('equal')
plt.plot(g.rbdry, g.zbdry, 'g--', label = 'Separatrix')
plt.plot(g.rlim, g.zlim, label = 'Limiter')
plt.title(r'Contour Plot of $B_R$')

plt.subplot(1,3,2)
plt.contourf(g.r_grid, g.z_grid, B_phi.T, cmap = 'plasma')
plt.colorbar(label = r'$B_{\phi}$')
plt.xlabel('R [m]')
plt.ylabel('Z [m]')
plt.axis('equal')
plt.plot(g.rbdry, g.zbdry, 'g--', label = 'Separatrix')
plt.plot(g.rlim, g.zlim, label = 'Limiter')
plt.title(r'Contour Plot of $B_{\phi}$')

plt.subplot(1,3,3)
plt.contourf(g.r_grid, g.z_grid, BZ, cmap = 'plasma')
plt.colorbar(label = r'$B_Z$')
plt.xlabel('R [m]')
plt.ylabel('Z [m]')
plt.axis('equal')
plt.plot(g.rbdry, g.zbdry, 'g--', label = 'Separatrix')
plt.plot(g.rlim, g.zlim, label = 'Limiter')
plt.title(r'Contour Plot of $B_Z$')

plt.legend()
plt.tight_layout()
plt.show()


# In[16]:


from scipy.interpolate import RegularGridInterpolator

BR_interp = RegularGridInterpolator((r_coords, z_coords), BR)
BZ_interp = RegularGridInterpolator((r_coords, z_coords), BZ)


# ## Using Euler

# In[17]:


r_start_points = np.linspace(g.rmaxis + 0.01, g.rbbbs.max()-0.01, 10)

r_trajectory = []
z_trajectory = []

z0 = g.zmaxis
h = 0.01

plt.figure(figsize=(6,6))

for i,initial_r_val in enumerate(r_start_points):

    current_r_val = initial_r_val
    current_z_val = z0

    r_trajectory.append([current_r_val])
    z_trajectory.append([current_z_val])

    revolution_completed = False


    prev_z_val = current_z_val

    while (not revolution_completed) and  \
          (r_coords.min() < current_r_val < r_coords.max()) and \
          (z_coords.min() < current_z_val < z_coords.max()):

        BR_local = BR_interp((current_r_val, current_z_val))
        BZ_local = BZ_interp((current_r_val, current_z_val))

        Bp_local = np.sqrt(BR_local**2 + BZ_local**2)

        if Bp_local == 0:
            break

        dR_dS = BR_local / Bp_local
        dZ_dS = BZ_local / Bp_local

        new_r_val = current_r_val + h * dR_dS
        new_z_val = current_z_val + h * dZ_dS


        if (prev_z_val < z0) and (new_z_val >= z0) :
            revolution_completed = True

        prev_z_val = current_z_val
        current_r_val = new_r_val
        current_z_val = new_z_val

        r_trajectory[i].append(current_r_val)
        z_trajectory[i].append(current_z_val)
        # step_count += 1

    plt.plot(r_trajectory[i], z_trajectory[i])

plt.plot(g.rbbbs, g.zbbbs, 'r--', alpha = 0.4, linewidth=2)
plt.plot(g.rmaxis, g.zmaxis, 'ro')
plt.plot(g.rlim, g.zlim, 'k', linewidth=2)
plt.axis('equal')
plt.xlabel("R")
plt.ylabel("Z")
plt.title("Poloidal Projection of Magnetic Field Lines")
plt.show()


# ## Using RK4

# In[18]:


def RK4(F, r_current, z_current, h):
  k1r, k1z = F(r_current, z_current)
  k2r, k2z = F(r_current + 0.5*h*k1r, z_current + 0.5*h*k1z)
  k3r, k3z = F(r_current + 0.5*h*k2r, z_current + 0.5*h*k2z)
  k4r, k4z = F(r_current + h*k3r, z_current + h*k3z)

  r_new = r_current + (h/6)*(k1r + 2*k2r + 2*k3r + k4r)
  z_new = z_current + (h/6)*(k1z + 2*k2z + 2*k3z + k4z)
  return r_new, z_new


# In[19]:


def F(r_current, z_current):
  BR_local = BR_interp((r_current, z_current))
  BZ_local = BZ_interp((r_current, z_current))
  BP = np.sqrt(BR_local**2 + BZ_local**2)
  if BP == 0.0:
    br = 0.0
    bz = 0.0
  else:
    br = BR_local/np.sqrt(BR_local**2 + BZ_local**2)
    bz = BZ_local/np.sqrt(BR_local**2 + BZ_local**2)
  return br, bz


# In[20]:


r_start = np.linspace(g.rmagx + 0.01, g.rbdry.max() - 0.01, 10)
z_start = g.zmaxis
h = 0.001

r_trajectory = []
z_trajectory = []

plt.figure(figsize=(6,6))

for i, r_initial in enumerate(r_start):
    r_current = r_initial
    z_current = z_start

    r_trajectory.append([r_current])
    z_trajectory.append([z_current])

    max_step = 50000
    step_count = 0

    revolution_completed = False
    while (not revolution_completed) and  \
          (r_coords.min() < current_r_val < r_coords.max()) and \
          (z_coords.min() < current_z_val < z_coords.max()) and \
          (step_count < max_step):

        step_count += 1

        r_new, z_new = RK4(F, r_current, z_current, h)

        if (z_current < z0) and (z_new >= z0) :
            revolution_completed = True

        r_current, z_current = r_new, z_new

        r_trajectory[i].append(r_new)
        z_trajectory[i].append(z_new)

    plt.plot(r_trajectory[i], z_trajectory[i], label=f"Line {i}")

plt.plot(g.rbbbs, g.zbbbs, 'r--', linewidth=2, alpha=0.4, label='Boundary')
plt.plot(g.rmaxis, g.zmaxis, 'ro', markersize=2, label='Magnetic Axis')
plt.plot(g.rlim, g.zlim, 'k', linewidth=2)
plt.axis('equal')
plt.xlabel("R")
plt.ylabel("Z")
plt.title("RK4 Poloidal Projection of Field Lines")
plt.legend()
plt.show()


# ## Using Constant Poloidal Flux for Adaptive RK4 tracing

# In[21]:


r_start = np.linspace(g.rmagx + 0.09, g.rbdry.max()-0.0001, 10)
z_start = g.zmaxis
tolerance = 1e-2
initial_h = -0.001

psi_interp = RegularGridInterpolator((r_coords, z_coords), g.psi)

r_trajectory = []
z_trajectory = []

plt.figure(figsize=(6,6))

for i, r_initial in enumerate(r_start):
    r_current = r_initial
    z_current = z_start

    psi_start = psi_interp(np.array([[r_current, z_current]]))[0]

    r_trajectory.append([r_current])
    z_trajectory.append([z_current])

    max_step = 50000
    step_count = 0


    while (step_count < max_step):

        psi_current = psi_interp(np.array([[r_current, z_current]]))[0]

        step_count += 1
        h_current_step = initial_h

        # Adaptive step refinement loop
        refined_step_found = False
        for _ in range(15):
            r_new_candidate, z_new_candidate = RK4(F, r_current, z_current, h_current_step)

            psi_new = psi_interp(np.array([[r_new_candidate, z_new_candidate]]))[0]

            if abs(psi_new - psi_start) < tolerance:
                r_new, z_new = r_new_candidate, z_new_candidate
                refined_step_found = True
                break
            else:
                h_current_step /= 2.0
                if h_current_step < 1e-7:
                    r_new, z_new = r_new_candidate, z_new_candidate
                    refined_step_found = True
                    break

        if step_count >1000:
          if  abs(r_initial - r_new) < 1e-4:
            break

        r_current, z_current = r_new, z_new

        r_trajectory[i].append(r_new)
        z_trajectory[i].append(z_new)



    plt.plot(r_trajectory[i], z_trajectory[i], label=f"Line {i}")

plt.plot(g.rbbbs, g.zbbbs, 'r--', linewidth=2, alpha=0.4, label='Boundary')
plt.plot(g.rmaxis, g.zmaxis, 'ro', markersize=2, label='Magnetic Axis')
plt.plot(g.rlim, g.zlim, 'k', linewidth=2)
plt.axis('equal')
plt.xlabel("R")
plt.ylabel("Z")
plt.title("Adaptive Poloidal Projection of Magnetic Field Lines")
plt.legend()
plt.show()


# In[22]:


len(r_trajectory[1])


# In[23]:


# r_grid = np.concatenate(r_trajectory)


# In[24]:


# r_grid.shape


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




