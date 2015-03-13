#Enrichment_Animation.py
#Use this module to get data from the Calculation and Animation modules for two different diffusivities, then plot them together (normalized to 1: being the starting concentration of the isotope) and animate them together. More importantly, this needs to subtract the two isotope concnetration gradients in order to define the enrichment along the (long) x-axis of the experiment. I would also like to put the "stop" at a point where the diffusion reaches the end of the grid, so It will automatically give the ned result that I need. 

import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
from Parameters import nx, nt, UL, UR, dx, UnL, UnR
from BC_Vector import Dsp, D_secsp, bc, bc_sec
from Initial_Condition import condition
import scipy as sp
from scipy.sparse import csr_matrix, linalg

def run_animate_two(timestep):
    global n
    global u
    global u_sec
    global nt
    if timestep > nt:
        return None
    U = u
    U_sec = u_sec
    U = np.delete(U, [0, nx-1], axis = 0)
    U_sec = np.delete(U_sec, [0, nx-1], axis = 0)
    U = U + bc
    U_sec = U_sec + bc_sec
    U = sp.sparse.linalg.spsolve(Dsp, U.T).T
    U_sec = sp.sparse.linalg.spsolve(D_secsp, U_sec.T).T
    u = np.append(np.append(U, UR)[::-1], UL)[::-1]
    u_sec = np.append(np.append(U_sec, UR)[::-1], UL)[::-1]
    enrich = abs(u - u_sec)
    line1.set_ydata(u)
    line2.set_ydata(u_sec)
    line3.set_ydata(enrich)
    return line1, line2, line3

x = condition[:, 1]
u = condition[:, 0]
u_sec = condition[:, 0]
en_data = np.empty(len(condition))
fig1 = plt.figure()
line1, line2, line3 = plt.plot(x, u, 'b-', x, u_sec, 'r-', x, en_data, 'y-')
plt.xlabel('Distance (x)')
plt.ylabel('Concentration (u)')
plt.title('Argon Isotope Concentrations')
ani = animation.FuncAnimation(fig1, run_animate_two)
plt.show()
exit()