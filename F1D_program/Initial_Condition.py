#Initial_Condition.py
from scipy import shape
from Grid import gridx
from Parameters import nx, dx, L
import numpy as np

#creates a 2d ?x2 matrix containing zeros in the first column and gridx in the second column. 
condition = np.zeros([nx, 2])
condition[:, 1] = np.array(gridx)

gridu = []
for i in gridx:
    if i >= L/2:
        u = 1.0
        #u = 0.0
    else:
        u = 0.0
    gridu.append(u)
condition[:, 0] = np.array(gridu)

