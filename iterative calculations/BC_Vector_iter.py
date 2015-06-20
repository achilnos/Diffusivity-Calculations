#BC_Vector_iter.py
#This module contains the script to create a vector "D" containing a Dirichlet boundary condition for this system. 

from Parameters_iter import nx, dt, dx, UL, UR, UnL, UnR, vis, vis_sec
import numpy as np
from scipy.sparse import csr_matrix

def BC_matrix_gen():
    vis, vis_sec = diffusion_magnitudes()    
    L, nx, nt, dt, dx, UL, UR, UnL, UnR = static_parameters()
    bc = np.zeros([nx-2, 1])
    bc_sec = np.zeros([nx-2, 1])
    #Dirichlet boundary conditions
    #bc[0] = vis*dt*UL/(dx**2)
    #bc_sec[0] = vis_sec*dt*UL/(dx**2)
    #bc[nx-3] = vis*dt*UR/(dx**2)
    #bc_sec[nx-3] = vis_sec*dt*UR/(dx**2)
    #bc = bc.T
    #bc_sec = bc_sec.T
    #Neumann boundary conditions
    bc[0] = vis*dt*-UnL/(dx) #N. B. C. left side
    bc_sec[0] = vis_sec*dt*-UnL/(dx) #N. B. C. left side
    bc[nx-3] = vis*dt*UnR/(dx) #N. B. C. right side
    bc_sec[nx-3] = vis_sec*dt*UnR/(dx)#N. B. C. right side
    bc = bc.T
    bc_sec = bc_sec.T
    #calculates the coefficent matrix for the implicent scheme. 
    A = np.delete(np.hstack((np.zeros([nx-2, 1]),np.identity(nx-2))), nx-2, 1)
    # Dirichlet B.C.s calcs
    #B = A + A.transpose() - 2 * np.identity(nx-2)
    #D = np.identity(nx-2) - (vis*dt/dx**2)*B
    #D_sec = np.identity(nx-2) - (vis_sec*dt/dx**2)*B
    #Dsp = csr_matrix(np.identity(nx-2) - (vis*dt/dx**2)*B)
    #D_secsp = csr_matrix(np.identity(nx-2) - (vis_sec*dt/dx**2)*B)
    # Neumann B. C. s calcs
    B = A + A.transpose() - 2 * np.identity(nx-2)
    B[0,0] = -1
    B[nx-3,nx-3] = -1
    D = np.identity(nx-2) - (vis*dt/dx**2)*B
    D_sec = np.identity(nx-2) - (vis_sec*dt/dx**2)*B
    Dsp = csr_matrix(np.identity(nx-2) - (vis*dt/dx**2)*B)
    D_secsp = csr_matrix(np.identity(nx-2) - (vis_sec*dt/dx**2)*B)

result = BC_matrix_gen()
print result#currently prints none!