#Plot.py
#plots a specific timestep from a text file
from scipy import shape
import numpy as np
from numpy import loadtxt
from Initial_Condition import condition
import matplotlib.pylab as plt

Arforty = np.loadtxt('Ar40','float128',)
Arthirtysix = np.loadtxt('Ar36','float128',)

def enrich_cor(Ar40, Ar36, x):
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

def Plot():    
    #endtime for T = 1100 is 13672
    a = Arforty[8000,:]
    b = Arthirtysix[8000,:]
    Ar40 = a[1:len(a)]
    Ar36 = b[1:len(b)]
    R = Ar36/Ar40
    c = condition[:, 1]
    x = c[1:len(b)]
    e, Caev, Cbev = enrich_cor(Ar40, Ar36, x)
    #plt.plot(x, Caev, 'b-', x, Cbev, 'r-')
    plt.plot(x, R)
    #plt.plot(x, e)
    #plt.plot(x, Caev)
    plt.xlabel('distance (m)')
    plt.ylabel('36/40 ratio (real, linear)')
    #plt.ylabel('Real concentration (Log)')
    #plt.ylabel('Enrichment 36Ar (permil)')
    #plt.ylabel('[36Ar] (ppm)')
    #plt.title('Ar isotope enrichment')
    plt.title('36Ar (blue) vs. 40Ar (red) in DC')
    #plt.title('36 Ar Concentration (Log Scale)')
    plt.grid(True)
    #plt.savefig("F1D_result.png")
    #plt.subplot(111).set_yscale('log')
    plt.show()
    exit()

Plot()
