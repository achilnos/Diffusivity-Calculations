#Grid.py
#This script defines the grid for F1D
#It imports "dx" from the F1D_Parameters module and uses a decimal range function to generate the grid and place this string into "grid". Formats the numbers to display as exponentials if very small. 

from Parameters import dx, L, nx
import numpy as np
from scipy import shape

def drange(start, stop, step):
     r = start
     while r <= stop:
     	yield r
     	r += step
     #yield r# not sure why I keep having to flip this on and off...
gridx = []
#gridx = np.empty([nx, 0])
["%g" % x for x in drange(0.0, L, dx)]
for x in drange(0.0, L, dx):
    gridx.append(x+L/(nx*2))
