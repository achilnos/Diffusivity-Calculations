#runtime_calc.py


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib as pp
#from Parameters import dx, L, nx
from scipy import shape
from scipy.sparse import csr_matrix
import scipy as sp
from scipy.sparse import csr_matrix, linalg
import csv
import math

#dynamic model parameters; creates the set of pressure and temperature values to use to calculate D. 

def temperature_set():
    temperature_low = 900
    temperature_high = 1200
    temperature_increment = 20
    temperature_vector = np.asarray(range(temperature_low, temperature_high, temperature_increment))
    return temperature_vector

def pressure_set():
    pressure_low = 50
    pressure_high = 100
    pressure_increment = 10
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
    nt = 1000000
    #Width of the time step (what are the units for this?)
    dt = 10
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
    return gridx
#sets up the boundary condition vector

def BC_matrix_gen(dif_40, dif_39, dif_38, dif_37, dif_36, beta_40, beta_39, beta_38, beta_37, beta_36):  
    #Common parameters
    L, nx, nt, dt, dx, UL, UR, UnL, UnR = static_parameters()
    bc_40 = np.zeros([nx-2, 1])
    bc_39 = np.zeros([nx-2, 1])
    bc_38 = np.zeros([nx-2, 1])
    bc_37 = np.zeros([nx-2, 1])
    bc_36 = np.zeros([nx-2, 1])    
    #Dirichlet boundary conditions
    bc_40[0] = dif_40*dt*UL/(dx**2)
    bc_39[0] = dif_39*dt*UL/(dx**2)
    bc_38[0] = dif_38*dt*UL/(dx**2)
    bc_37[0] = dif_37*dt*UL/(dx**2)
    bc_36[0] = dif_36*dt*UL/(dx**2)
    bc_40[nx-3] = dif_40*dt*UR/(dx**2)
    bc_39[nx-3] = dif_39*dt*UR/(dx**2)
    bc_38[nx-3] = dif_38*dt*UR/(dx**2)
    bc_37[nx-3] = dif_37*dt*UR/(dx**2)
    bc_36[nx-3] = dif_36*dt*UR/(dx**2)
    #Neumann boundary conditions
    #bc_40[0] = dif_40*dt*-UnL/(dx) #sets N. B. C. left side (40Ar)
    #bc_39[0] = dif_39*dt*-UnL/(dx) #sets N. B. C. left side (39Ar)
    #bc_38[0] = dif_38*dt*-UnL/(dx) #sets N. B. C. left side (38Ar)
    #bc_37[0] = dif_37*dt*-UnL/(dx) #sets N. B. C. left side (37Ar)
    #bc_36[0] = dif_36*dt*-UnL/(dx) #sets N. B. C. left side (36Ar)
    #bc_40[nx-3] = dif_40*dt*UnR/(dx) #sets N. B. C. right side (40Ar)
    #bc_39[nx-3] = dif_39*dt*UnR/(dx)#sets N. B. C. right side (39Ar)
    #bc_38[nx-3] = dif_38*dt*UnR/(dx)#sets N. B. C. right side (39Ar)
    #bc_37[nx-3] = dif_37*dt*UnR/(dx)#sets N. B. C. right side (39Ar)
    #bc_36[nx-3] = dif_36*dt*UnR/(dx)#sets N. B. C. right side (39Ar)    
    #matrix operation for Neumann or Dirichlet B.C.s
    bc_40 = bc_40.T
    bc_39 = bc_39.T
    bc_38 = bc_38.T
    bc_37 = bc_37.T
    bc_36 = bc_36.T
    #calculates the coefficent matrix for the implicent scheme. 
    A = np.delete(np.hstack((np.zeros([nx-2, 1]),np.identity(nx-2))), nx-2, 1)
    B = A + A.transpose() - 2 * np.identity(nx-2)    
    #matrix operations for Neumann or Dirichlet B.C.s
    B = A + A.transpose() - 2 * np.identity(nx-2)
    D_40 = np.identity(nx-2) - (dif_40*dt/dx**2)*B
    D_39 = np.identity(nx-2) - (dif_39*dt/dx**2)*B
    D_38 = np.identity(nx-2) - (dif_38*dt/dx**2)*B
    D_37 = np.identity(nx-2) - (dif_37*dt/dx**2)*B
    D_36 = np.identity(nx-2) - (dif_36*dt/dx**2)*B
    B = A + A.transpose() - 2 * np.identity(nx-2)    
    #Neumann only B. C. s calcs
    #B[0,0] = -1
    #B[nx-3,nx-3] = -1
    #matrix operations for Neumann or Dirichlet B.C.s
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

#sets the intial condition (currently set to the empty condition, rather than the step function)

def initial_condition():
    gridx = grid_maker()
    L, nx, nt, dt, dx, UL, UR, UnL, UnR = static_parameters()
    condition = np.zeros([nx, 2])
    condition[:, 1] = np.array(gridx)
    gridu = []
    for i in gridx:
        if i >= L/2:
            #u = 1.0
            u = 0.0
        else:
            u = 0.0
        gridu.append(u)
    condition[:, 0] = np.array(gridu)
    return condition

#calculation

def Calculator(dif_40, dif_39, dif_38, dif_37, dif_36, beta_40, beta_39, beta_38, beta_37, beta_36):#what is beta for? I should remove it...
    L, nx, nt, dt, dx, UL, UR, UnL, UnR = static_parameters()
    D_40_sp, D_39_sp, D_38_sp, D_37_sp, D_36_sp, bc_40, bc_39, bc_38, bc_37, bc_36 = BC_matrix_gen(dif_40, dif_39, dif_38, dif_37, dif_36, beta_40, beta_39, beta_38, beta_37, beta_36)
    Dsp = D_40_sp#temporary until I add the rest of the diffusivities to the calculation
    bc = bc_40#temporary until I add the rest of the bc vectors to the calculation
    condition = initial_condition()
    u = (condition[:, 0])
    for timestep in range(0, nt):
        #print timestep
        if u[88]*10. > u[nx-1] :
            #print 'Runtime to equilibrium is ', timestep*dt, ' seconds ', 'or ', timestep*dt/60./60., ' hours.'
        #if u[1] > 0.00001:
        #    print 'Optimal runtime is', timestep*dt, 'seconds.'
            return (u, timestep * dt)
        U = u
        U = np.delete(U, [0, nx-1], axis = 0)
        U = U + bc
        U = sp.sparse.linalg.spsolve(Dsp, U.T).T
        u = np.append(np.append(U, UR)[::-1], UL)[::-1]
        #u = np.append(np.append(U, U[0]-UnL*dx)[::-1], U[nx-3]+UnR*dx)[::-1]
    print 'alloted computation time exceeded'
    return (u, nt * dt)

#iterator and plotting function
#for the given input of pressure and temperature, finds a value of viscosity

def diffusion_magnitudes():
    with open('200_10.txt', 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=' ')
        L, nx, nt, dt, dx, UL, UR, UnL, UnR = static_parameters()
        temperatures, pressures = conditions_set()
        #Water Concentration (weight percent)
        W = 0.04#RLS-41
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
            u, time_to_stop = Calculator(dif_40, dif_39, dif_38, dif_37, dif_36, beta_40, beta_39, beta_38, beta_37, beta_36)
            T = str(int(T))
            P = str(int(P))
            dif_pri = str(dif_40)
            #dif_sec = str(dif_39)
            #dif_tert = str(dif_38)
            #dif_quat = str(dif_37)
            #dif_sept = str(dif_36)
            seconds = str(time_to_stop)
            minutes = str(time_to_stop/60.)
            hours = str(time_to_stop/60./60.)
            hours_time = str(int(math.floor(float(hours))))
            remainder = (float(hours) - math.floor(float(hours)))
            minutes_time = str(int(math.floor(float(remainder)*60.)))
            remainder = (float(minutes_time) - math.floor(float(minutes_time)))
            seconds_time = str(int(math.floor(float(remainder)*60.)))
            runtime = hours_time, minutes_time, seconds_time
            #pads with left zero if single digit to get HH:MM:SS format
            i = 0
            runtime = list(runtime)
            while i < 3:
                n = str(runtime[i])
                if len(n) < 2:
                    runtime[i] = n.zfill(2)
                i = i + 1
            runtime = tuple(runtime)
            runtime = ":".join(runtime)
            output = dif_pri, T, P, seconds, minutes, hours, runtime 
            output = ','.join(output)
            print output
            spamwriter.writerow(output)

#plotting step here
diffusion_magnitudes()
