#Animation.py
#This script calculates the result of each timestep and plots the result at the final timestep. An optional animation for eachtimestep is included and can be used if the lines () are commented in.  

import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
from F1D_Parameters import nx, nt, UL, UR, dx, UnL, UnR
from F1D_BC_Vector import Dsp, bc
from F1D_Initial_Condition import condition
import scipy as sp
from scipy.sparse import csr_matrix, linalg

def run_animate(timestep):
    global u
    global nt
    if timestep > nt:
        return None
    else: 
        print timestep
    U = u
    U = np.delete(U, [0, nx-1], axis = 0)
    U = U + bc
    U = sp.sparse.linalg.spsolve(Dsp, U.T).T
    u = np.append(np.append(U, UR)[::-1], UL)[::-1]
    graph.set_ydata(u)
    return graph,


x = condition[:, 1]
u = condition[:, 0]
fig1 = plt.figure()
graph, = plt.plot(x, u)
plt.xlabel('Distance (x)')
plt.ylabel('Concentration (u)')
plt.title('Argon Concentration')
ani = animation.FuncAnimation(fig1, run_animate)
plt.show()
exit()


