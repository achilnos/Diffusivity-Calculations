#Enrichment_Calc.py

import numpy as np
import matplotlib.pylab as plt
import matplotlib.animation as animation
from Parameters import nx, nt, UL, UR, dx, UnL, UnR, vis, vis_sec
from BC_Vector import Dsp, D_secsp, bc, bc_sec
from Initial_Condition import condition
import scipy as sp
from scipy.sparse import csr_matrix, linalg

def Calculator(vis, vis_sec):
    u = (condition[:, 0])
    u_sec = condition[:, 0]
    for timestep in range(0, nt):
        if u[1] > 0.00001:
            print 'Optimal runtime is', timestep, 'seconds.'
            return u, u_sec
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
    return u, u_sec

def enrich_cor(vis, vis_sec):
    Ar36, Ar40 = Calculator(vis, vis_sec)
    #afrac is the initial fraction of 36Ar in the argon enriched couple at the begining of the experiment. Currently using atmospheric argon values. 
    afrac = 0.00337
    #bfrac is the initial fraction of 40Ar in the argon enriched couple at the begining of the experiment. Currently using atmospheric argon values. 
    bfrac = 0.996
    #Cin is the intial cocentration of argon per unit length in the enriched diffusion couple. In this formulation this is in parts per million. Using 350 ppm
    Cin = 350.0
    #Rin is the ratio of 36 to 40 Ar (initial)
    Rin = afrac/bfrac
    #Cain is the initial concentration of argon 36 in the DC (ppm). 
    Cain = Cin*afrac
    #Cbin is the initial concentration of argon 40 in the DC (ppm). 
    Cbin = Cin*bfrac
    #Faev is the evolved fraction of argon 36 at some point in x and t during the experiment.
    Faev = Ar36*afrac
    #Fbev is the evolved fraction of argon 40 at some point in x and t during the experiment.
    Fbev = Ar40*bfrac
    #Caev is the evolved concentration of 36Ar along the x axis at a particular timestep.
    Caev = Ar36*Cain
    #Cbev is the evolved concentration of 40Ar along the x axis at a particular timestep
    Cbev = Ar40*Cbin
    #Rev is the evolved ratio of 36Ar to 40Ar along the x axis at a particular timestep.    
    Caev = Caev[1:len(Ar40)]
    Cbev = Cbev[1:len(Ar36)]
    Rev = Caev/Cbev#invalid value encountered in divide
    #e is the enrichment (del) relative to atmospheric in per mil. 
    permil = (Rev/(Rin-1))/(Rin*1000.0)
    return permil, Caev, Cbev

def Plot(vis, vis_sec):    
    permil, Caev, Cbev = enrich_cor(vis, vis_sec)
    x = condition[:, 1]
    x = x[1:len(x)]
    f1 = plt.figure()
    plt.xlabel('distance (m)')
    plt.ylabel('Concentration (ppm)')
    plt.title('Argon Isotope Concentrations')
    plt.subplot(111).set_yscale('log')
    plt.grid(True)
    f2 = plt.figure()
    ax1 = f1.add_subplot(111)
    ax2 = f2.add_subplot(111)
    ax1.plot(x, Caev, 'r-', x, Cbev, 'b-')
    ax2.plot(x, permil, 'b-')
    plt.xlabel('distance (m)')
    plt.ylabel('Isotope enrichment (permil)')
    plt.title('Argon Isotope Enrichment')
    plt.grid(True)
    
    plt.show()

Plot(vis,vis_sec)