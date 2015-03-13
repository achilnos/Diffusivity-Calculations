#Plot_iter.py
#the purpose of this program is to generate a matrix of the equilibrium runtimes and plot them as a as a color contour plot, T = x-axis, P = y-axis, t = color gradient. 

from Calculation_iter import times_to_equilibrium
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from Parameters_iter import P_T
P_T_t = np.vstack((P_T,times_to_equilibrium))


def getMeshGrid(P_T_t):#not currently utilized. used to process data for surface plot thing. 
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
X = P_T_t[:,1]
Y = P_T_t[:,2]
Z = P_T_t[:,3]

np.savetxt('runtimePT.csv', P_T_t, delimiter=',')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(X,Y,Z)
ax.plot_surface(P_T_t[:,0], P_T_t[:,1], P_T_t[:,2])
ax.plot_trisurf(P_T_t[:,0], P_T_t[:,1], P_T_t[:,2])#surf = ax.plot_surface( X, Y, Z, rstride=1, cstride=1, #cmap=cm.coolwarm,
plt.show()
#        linewidth=0, antialiased=False)

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

#fig.colorbar(surf, shrink=0.5, aspect=5)
#plt.show()