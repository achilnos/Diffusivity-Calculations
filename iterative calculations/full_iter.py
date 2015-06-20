#full_iter.py


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib as pp
#from Parameters import dx, L, nx
from scipy import shape
from scipy.sparse import csr_matrix

#dynamic model parameters; creates the set of pressure and temperature values to use to calculate D. 

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
    return np.reshape(temperature_array, len(temperature_vector) * len(pressure_vector)), np.reshape(pressure_array, len(temperature_vector) * len(pressure_vector))

#static (not iterating) model parameters

def static_parameters():
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
    #Left Dirichlet B.C. (y'' + y = 0)
    UL = 0
    #Right Boundary condition is Dirichlet B.C: (y'' + y = 0)
    UR = 1
    #Left Neumann B.C.
    UnL = 0
    #Left Neumann B.C.
    UnR = 1
    return L, nx, nt, dt, dx, UL, UR, UnL, UnR

#for the given input of pressure and temperature, finds a value of viscosity

def diffusion_magnitudes():
    L, nx, nt, dt, dx, UL, UR, UnL, UnR = static_parameters()
    temperatures, pressures = conditions_set()
    #Water Concentration (weight percent)
    W = 0.1669364#not sure if this is the correct value!!!
    #mass slection coefficent
    E = 0.43
    for i in range(len(temperatures)):
        #Temperature (Kelvin)
        T = temperatures[i]
        #Pressure (MPa)
        P = pressures[i]
        #diffusivities for 40Ar, 39Ar, 38Ar, 37Ar, 36Ar
        dif_40 = (np.exp((14.627-17913/T-2.569*P/T)+(35936/T+27.42*P/T)*W))*10**-12
        dif_39 = dif_40*(40.0/39.0)**(E/2.0)
        dif_38 = dif_40*(40.0/38.0)**(E/2.0)
        dif_37 = dif_40*(40.0/37.0)**(E/2.0)
        dif_36 = dif_40*(40.0/36.0)**(E/2.0)
        #stability criterion beta ...what is the stability criterion anyway?
        beta_40 = dif_40*dt/(dx*dx)
        beta_39 = dif_39*dt/(dx*dx)
        beta_38 = dif_38*dt/(dx*dx)
        beta_37 = dif_37*dt/(dx*dx)
        beta_36 = dif_36*dt/(dx*dx)
    return dif_40, dif_39, dif_38, dif_37, dif_36, beta_40, beta_39, beta_38, beta_37, beta_36

#sets up a mesh for the model

def drange(start, stop, step):
     r = start
     while r <= stop:
     	yield r
     	r += step#yield r# not sure why I keep having to flip this on and off...
     
def grid_maker():
    L, nx, nt, dt, dx, UL, UR, UnL, UnR = static_parameters()
    gridx = []
    #gridx = np.empty([nx, 0])
    ["%g" % x for x in drange(0.0, L, dx)]
    for x in drange(0.0, L, dx):
        gridx.append(x+L/(nx*2))   

#sets up the boundary condition vector

def BC_matrix_gen():
    dif_40, dif_39, dif_38, dif_37, dif_36, beta_40, beta_39, beta_38, beta_37, beta_36 = diffusion_magnitudes()    
    vis = dif_40
    vis_sec = diff_36
    L, nx, nt, dt, dx, UL, UR, UnL, UnR = static_parameters()
    bc_40 = np.zeros([nx-2, 1])
    bc_39 = np.zeros([nx-2, 1])
    bc_38 = np.zeros([nx-2, 1])
    bc_37 = np.zeros([nx-2, 1])
    bc_36 = np.zeros([nx-2, 1])
    #Dirichlet boundary conditions
    #bc_40[0] = dif_40*dt*UL/(dx**2)
    #bc_39[0] = dif_39*dt*UL/(dx**2)
    #bc_38[0] = dif_38*dt*UL/(dx**2)
    #bc_37[0] = dif_37*dt*UL/(dx**2)
    #bc_36[0] = dif_36*dt*UL/(dx**2)
    #bc_40[nx-3] = dif_40*dt*UR/(dx**2)
    #bc_39[nx-3] = dif_39*dt*UR/(dx**2)
    #bc_38[nx-3] = dif_38*dt*UR/(dx**2)
    #bc_37[nx-3] = dif_37*dt*UR/(dx**2)
    #bc_36[nx-3] = dif_36*dt*UR/(dx**2)
    #bc_40 = bc_40.T
    #bc_39 = bc_39.T
    #bc_38 = bc_38.T
    #bc_37 = bc_37.T
    #bc_36 = bc_36.T
    #Neumann boundary conditions
    bc_40[0] = dif_40*dt*-UnL/(dx) #sets N. B. C. left side (40Ar)
    bc_39[0] = dif_39*dt*-UnL/(dx) #sets N. B. C. left side (39Ar)
    bc_38[0] = dif_38*dt*-UnL/(dx) #sets N. B. C. left side (38Ar)
    bc_37[0] = dif_37*dt*-UnL/(dx) #sets N. B. C. left side (37Ar)
    bc_36[0] = dif_36*dt*-UnL/(dx) #sets N. B. C. left side (36Ar)
    bc_40[nx-3] = dif_40*dt*UnR/(dx) #sets N. B. C. right side (40Ar)
    bc_39[nx-3] = dif_39*dt*UnR/(dx)#sets N. B. C. right side (39Ar)
    bc_38[nx-3] = dif_38*dt*UnR/(dx)#sets N. B. C. right side (39Ar)
    bc_37[nx-3] = dif_37*dt*UnR/(dx)#sets N. B. C. right side (39Ar)
    bc_36[nx-3] = dif_36*dt*UnR/(dx)#sets N. B. C. right side (39Ar)
    bc_40 = bc_40.T
    bc_39 = bc_39.T
    bc_38 = bc_38.T
    bc_37 = bc_37.T
    bc_36 = bc_36.T
    #calculates the coefficent matrix for the implicent scheme. 
    A = np.delete(np.hstack((np.zeros([nx-2, 1]),np.identity(nx-2))), nx-2, 1)
    # Dirichlet B.C.s calcs
    #B = A + A.transpose() - 2 * np.identity(nx-2)
    #D_40 = np.identity(nx-2) - (dif_40*dt/dx**2)*B
    #D_39 = np.identity(nx-2) - (dif_39*dt/dx**2)*B
    #D_38 = np.identity(nx-2) - (dif_38*dt/dx**2)*B
    #D_37 = np.identity(nx-2) - (dif_37*dt/dx**2)*B
    #D_36 = np.identity(nx-2) - (dif_36*dt/dx**2)*B
    #D_40_sp = csr_matrix(np.identity(nx-2) - (dif_40*dt/dx**2)*B)
    #D_39_sp = csr_matrix(np.identity(nx-2) - (dif_39*dt/dx**2)*B)
    #D_38_sp = csr_matrix(np.identity(nx-2) - (dif_38*dt/dx**2)*B)
    #D_37_sp = csr_matrix(np.identity(nx-2) - (dif_37*dt/dx**2)*B)
    #D_36_sp = csr_matrix(np.identity(nx-2) - (dif_36*dt/dx**2)*B)
    # Neumann B. C. s calcs
    B = A + A.transpose() - 2 * np.identity(nx-2)
    B[0,0] = -1
    B[nx-3,nx-3] = -1
    D_40 = np.identity(nx-2) - (dif_40*dt/dx**2)*B
    D_39 = np.identity(nx-2) - (dif_39*dt/dx**2)*B
    D_38 = np.identity(nx-2) - (dif_38*dt/dx**2)*B
    D_37 = np.identity(nx-2) - (dif_37*dt/dx**2)*B
    D_36 = np.identity(nx-2) - (dif_36*dt/dx**2)*B
    D_40_sp = csr_matrix(np.identity(nx-2) - (dif_40*dt/dx**2)*B)
    D_39_sp = csr_matrix(np.identity(nx-2) - (dif_39*dt/dx**2)*B)
    D_38_sp = csr_matrix(np.identity(nx-2) - (dif_38*dt/dx**2)*B)
    D_37_sp = csr_matrix(np.identity(nx-2) - (dif_37*dt/dx**2)*B)
    D_36_sp = csr_matrix(np.identity(nx-2) - (dif_36*dt/dx**2)*B)
    return D_40_sp, D_39_sp, D_38_sp, D_37_sp, D_36_sp, bc_40, bc_39, bc_38, bc_37, bc_36

