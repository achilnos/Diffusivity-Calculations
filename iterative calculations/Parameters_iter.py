#Parameters_iter.py

import numpy as np

def KIF_Mag():
    #Temperature (Kelvin)
    T = 1000
    #Pressure (MPa)
    P = 100
    #Water Concentration (weight percent)
    #W = 0.1669364
    W = 0.01
    #mass slection coefficent
    E = 0.43
    vis = (np.exp((14.627-17913/T-2.569*P/T)+(35936/T+27.42*P/T)*W))*10**-12
    vis_sec = vis*(40.0/36.0)**(E/2.0)
    return vis, vis_sec

#Length of System (meters)
L = 0.00169
#Number of steps in space
nx = 100
#Number of steps in time
nt = 9999999
#Width of the time step (what are the units for this?)
dt = 1
#Width of the space step
dx = L / ( nx - 1.0 )
#diffusivity
vis, vis_sec = KIF_Mag()
#stability criterion
beta = vis*dt/(dx*dx)
beta_sec = vis_sec*dt/(dx*dx)
#Left Dirichlet B.C. (y'' + y = 0)
UL = 0
#Right Boundary condition is Dirichlet B.C: (y'' + y = 0)
UR = 1
#Left Neumann B.C.
UnL = 0
#Left Neumann B.C.
UnR = 1

print "parameter accessed"
