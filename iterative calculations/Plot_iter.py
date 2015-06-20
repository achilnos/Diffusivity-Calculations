#Plot_iter.py
#the purpose of this program is to generate a matrix of the equilibrium runtimes and plot them as a as a color contour plot, T = x-axis, P = y-axis, t = color gradient. 

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib as pp

def temperature_set():
    temperature_low = 900
    temperature_high = 1100
    temperature_increment = 25
    temperature_vector = np.asarray(range(temperature_low, temperature_high, temperature_increment))
    return temperature_vector

def pressure_set():
    pressure_low = 50
    pressure_high = 100
    pressure_increment = 5
    pressure_vector = np.asarray(range(pressure_low, pressure_high, pressure_increment))
    return pressure_vector

def conditions_set():
    temperature_vector = temperature_set()
    pressure_vector = pressure_set()
    temperature_array = np.empty([len(temperature_vector), len(pressure_vector)])
    pressure_array = np.empty([len(temperature_vector), len(pressure_vector)])
    for i in range(len(temperature_vector)):
        temperature_array[i,:] = temperature_vector[i]
    for j in range(len(pressure_vector)):
        pressure_array[:,j] = pressure_vector[j]
    return temperature_array, pressure_array




#def getMeshGrid(P_T_t):#not currently utilized. used to process data for surface plot thing. Actually yes use this it worked. 
#    dataArray = P_T_t.T
#    """gets 2d coordinate grid and Z values in meshgrid format. requires values in
#    dataArray to have a rectangular region of x-y space covered uniformly"""
#    xs = dataArray[:,0]
#    ys = dataArray[:,1]
#    xmin,xmax = xs.min(), xs.max()
#    xstep = xs[xs!=xmin].min()-xmin
#    ymin,ymax = ys.min(), ys.max()
#    ystep = ys[ys!=ymin].min()-ymin
#    X = np.arange(xmin, xmax+xstep, xstep)
#    Y = np.arange(ymin, ymax+ystep, ystep)
#    X,Y = np.meshgrid(X,Y)
#    Z = np.zeros(X.shape)
#    height, width = X.shape
#    for i in range(0, height):
#        for j in range(0,width):
#            halfway = dataArray[dataArray[:,0]==X[i,j]] # finds all with that value of x
#            row = halfway[halfway[:,1]==Y[i,j]] # finds y value 
#            Z[i,j] = row[0,2]
#    return X,Y,Z

#P_T_t = P_T_t.T
#np.savetxt('runtimePT.txt', P_T_t, delimiter=',')
#X = P_T_t[:,0]
#Y = P_T_t[:,1]
#Z = P_T_t[:,2]
#print "length X is ",len(X), "length Y is ",len(Y), "length Z is ",len(Z)#this isn't the complete dataset!
#print "shape P_T_t is", np.shape(P_T_t)
#
#try using mayavi instead. It might be simpler to use than matplotlib for 3d. 
#
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#X = np.arange(-5, 5, 0.25)
#Y = np.arange(-5, 5, 0.25)
#X, Y = np.meshgrid(X, Y)#won't work with three vars. need to use a sci py interpolation function to get the relationship with Z. 
#scatty = ax.scatter(X, Y, Z)
#ax.set_zlim(-1.01, 1.01)
#
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#
#fig.colorbar(scatty, shrink=0.5, aspect=5)
#
#P = pp.imshow(Matrix, origin='lower', extent=(0, half_span * 2, 0, half_span * 2), cmap='YlGnBu', aspect='auto', interpolation='nearest')
#
#
#plt.show()