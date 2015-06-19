#Plot_iter.py
#the purpose of this program is to generate a matrix of the equilibrium runtimes and plot them as a as a color contour plot, T = x-axis, P = y-axis, t = color gradient. 

from Calculation_iter import times_to_equilibrium
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from Parameters_iter import P_T
import matplotlib as pp
P_T_t = np.vstack((P_T,times_to_equilibrium))

#################

import numpy as np

def KIF_Mag():
    #Temperature (Kelvin)
    #T = 1000
    T_low = 1000
    T_high = 1300
    T_incr = 3
    T_range = range(T_low, T_high, T_incr)
    #Pressure (MPa)
    #P = 200
    P_low = 170
    P_high = 270
    P_incr = 1
    P_range = range(P_low, P_high, T_incr)
    #Water Concentration (weight percent)
    W = 0.1669364#this may not be accurate. 
    #mass slection coefficent
    E = 0.43#this may not be accurate. 
#the mistake is here. sorted wrong. Combine P_T into a tuple to get them through this dude. Only one loop. Actually maybe not. ask Simon. 
    ########################
    P_range = list(range(P_low, P_high, T_incr))
    T_range = list(range(T_low, P_high, P_high))
    viscosity = np.empty((len(T_range), len(P_range)), dtype=float)
    for i in len(P_range):
        for j in len(T_range):
            vis =  (np.exp((14.627-17913/T-2.569*P/T)+(35936/T+27.42*P/T)*W))*10**-12
            viscosity[i,j] = vis
    return viscosity, P_range, T_range
    
    ########################
    vis_list = np.empty([2,len(T_range)*len(P_range)], dtype = float)
    P_T = np.empty([2,len(T_range)*len(P_range)], dtype = float)
    count = 0
    print P_range, T_range
    assert False
    for P in P_range:#this loop wrong?
        for T in T_range:
            vis = (np.exp((14.627-17913/T-2.569*P/T)+(35936/T+27.42*P/T)*W))*10**-12
            vis_sec = vis*(40.0/36.0)**(E/2.0)
            vis_list[0,count] = vis
            vis_list[1,count] = vis_sec
            P_T[0,count] = P
            P_T[1,count] = T
            count = count + 1#put the count in the loop?
    return vis_list, P_T

D, P_T = KIF_Mag()


#################

def getMeshGrid(P_T_t):#not currently utilized. used to process data for surface plot thing. Actually yes use this it worked. 
    dataArray = P_T_t.T
    """gets 2d coordinate grid and Z values in meshgrid format. requires values in
    dataArray to have a rectangular region of x-y space covered uniformly"""
    xs = dataArray[:,0]
    ys = dataArray[:,1]
    xmin,xmax = xs.min(), xs.max()
    xstep = xs[xs!=xmin].min()-xmin
    ymin,ymax = ys.min(), ys.max()
    ystep = ys[ys!=ymin].min()-ymin
    X = np.arange(xmin, xmax+xstep, xstep)
    Y = np.arange(ymin, ymax+ystep, ystep)
    X,Y = np.meshgrid(X,Y)
    Z = np.zeros(X.shape)
    height, width = X.shape
    for i in range(0, height):
        for j in range(0,width):
            halfway = dataArray[dataArray[:,0]==X[i,j]] # finds all with that value of x
            row = halfway[halfway[:,1]==Y[i,j]] # finds y value 
            Z[i,j] = row[0,2]
    return X,Y,Z

P_T_t = P_T_t.T
np.savetxt('runtimePT.txt', P_T_t, delimiter=',')
X = P_T_t[:,0]
Y = P_T_t[:,1]
Z = P_T_t[:,2]
print "length X is ",len(X), "length Y is ",len(Y), "length Z is ",len(Z)#this isn't the complete dataset!
print "shape P_T_t is", np.shape(P_T_t)

#try using mayavi instead. It might be simpler to use than matplotlib for 3d. 

fig = plt.figure()
ax = fig.gca(projection='3d')
#X = np.arange(-5, 5, 0.25)
#Y = np.arange(-5, 5, 0.25)
X, Y = np.meshgrid(X, Y)#won't work with three vars. need to use a sci py interpolation function to get the relationship with Z. 
scatty = ax.scatter(X, Y, Z)
#ax.set_zlim(-1.01, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(scatty, shrink=0.5, aspect=5)

P = pp.imshow(Matrix, origin='lower', extent=(0, half_span * 2, 0, half_span * 2), cmap='YlGnBu', aspect='auto', interpolation='nearest')


#plt.show()