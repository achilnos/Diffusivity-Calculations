#Parameters_iter.py
#Use this module to set the parameters for the model. 
#This script hold the parameters for a method to solve 1D Fouriers equation using an upwind scheme with an infite step function centered at zero as the intial condition. 
#Dir should be /Users/Nicholas/Documents/Python/Argon
#An example of acceptable experimental parameters: L = 0.008 m, nt = 36000 s (10h)
#T = 600C --> vis = 2.66384*10**-11
#T = 700C --> vis = 8.54662*10**-11
#T = 800C --> vis = 2.20661*10**-10
#T = 900C --> vis = 4.84651*10**-10
#T = 1000C --> vis = 2.39408210334*10**-10
#If computation can't do it try at 10 or 100 times speed by rasing vis and lowing nt. Currently at 100 times
#has been modified to reflect the radius of the current TZM sample. Currently using this to estimate the time until equilbration of the cylinder. 

import numpy as np

def KIF_Mag():
    #Temperature (Kelvin)
    #T = 1000
    T_low = 1200
    T_high = 1300
    T_incr = 5
    T_range = range(T_low, T_high, T_incr)
    #Pressure (MPa)
    #P = 200
    P_low = 800
    P_high = 900
    P_incr = 5
    P_range = range(P_low, P_high, T_incr)
    #Water Concentration (weight percent)
    W = 0.1669364#this may not be accurate. 
    #mass slection coefficent
    E = 0.43#this may not be accurate. 
#the mistake is here. sorted wrong. Combine P_T into a tuple to get them through this dude. Only one loop. Actually maybe not. ask Simon. 
    vis_list = np.empty([2,len(T_range)*len(P_range)], dtype = float)
    P_T = np.empty([2,len(T_range)*len(P_range)], dtype = float)
    count = 0
    for P in P_range:
        for T in T_range:
            vis = (np.exp((14.627-17913/T-2.569*P/T)+(35936/T+27.42*P/T)*W))*10**-12
            vis_sec = vis*(40.0/36.0)**(E/2.0)
            vis_list[0,count] = vis
            vis_list[1,count] = vis_sec
            P_T[0,count] = P
            P_T[1,count] = T
            count = count + 1
    return vis_list, P_T

D, P_T = KIF_Mag()
#Length of System
L = 0.00169
#Number of steps in space
nx = 1000
#Number of steps in time
nt = 360000
#Width of the time step (what are the units for this?)
dt = 1
#Width of the space step
dx = L / ( nx - 1.0 )
#Diffusion coefficent
#vis = 2.66384*10**-9
#vis_sec = 2.0*10**-9
#stability criterion
#beta = vis*dt/(dx*dx)#what are these for?
#beta_sec = vis*dt/(dx*dx)#what are these for?
#Left Dirichlet B.C. (y'' + y = 0)
UL = 0
#Right Boundary condition is Dirichlet B.C: (y'' + y = 0)
UR = 1
#Left Neumann B.C.
UnL = 0
#Left Neumann B.C.
UnR = 1

print "parameters generated"
