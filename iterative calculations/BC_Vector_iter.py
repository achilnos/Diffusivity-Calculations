#BC_Vector_iter.py
#This module contains the script to create a vector "D" containing a Dirichlet boundary condition for this system. 

from Parameters_iter import nx, dt, dx, UL, UR, UnL, UnR
import numpy as np
from scipy.sparse import csr_matrix


def BC_matrix_gen(nx, dt, dx, UL, UR, UnL, UnR):
    
    a, b = np.split(D, 2)#spit D into a vis and avis_sec array
    vis_values = a[0]
    vis_sec_values = b[0]
    #interate of vis and vis_sec arrays
#bc is a coefficent matrix for the implict scheme
    bc = np.zeros([nx-2, 1])
    bc_sec = np.zeros([nx-2, 1])
    Dsp_matrix = []
    D_secsp_matrix = []
    bc_matrix = []
    #bc_matrix = np.empty([len(vis_values), 1])
    bc_sec_matrix = []
    #bc_sec_matrix = np.empty([len(vis_values), 1])
# Neumann boundary conditions
    counter = 0
    for vis in vis_values:  
        bc[0] = vis*dt*-UnL/(dx)
        bc[nx-3] = vis*dt*UnR/(dx) 
        bc = bc.T
        A = np.delete(np.hstack((np.zeros([nx-2, 1]),np.identity(nx-2))), nx-2, 1)
        B = A + A.transpose() - 2 * np.identity(nx-2)
        B[0,0] = -1
        B[nx-3,nx-3] = -1
        D = np.identity(nx-2) - (vis*dt/dx**2)*B
        Dsp = csr_matrix(np.identity(nx-2) - (vis*dt/dx**2)*B)
        Dsp_matrix.append(Dsp) 
        bc_matrix.append(bc)   
        bc = bc.T
        counter = counter + 1
    counter = 0
    for vis_sec in vis_sec_values:
        bc_sec[0] = vis_sec*dt*-UnL/(dx) #N. B. C. left side
        bc_sec[nx-3] = vis_sec*dt*UnR/(dx)#N. B. C. right side
        bc_sec = bc_sec.T
        A = np.delete(np.hstack((np.zeros([nx-2, 1]),np.identity(nx-2))), nx-2, 1)
        B = A + A.transpose() - 2 * np.identity(nx-2)
        B[0,0] = -1
        B[nx-3,nx-3] = -1
        D_sec = np.identity(nx-2) - (vis_sec*dt/dx**2)*B
        D_secsp = csr_matrix(np.identity(nx-2) - (vis_sec*dt/dx**2)*B)
        D_secsp_matrix.append(D_secsp)
        bc_sec_matrix.append(bc_sec)
        bc_sec = bc_sec.T
    result = (bc_matrix, bc_sec_matrix, Dsp_matrix, D_secsp_matrix)
    return result
print 'boundary condition generated'
result = BC_matrix_gen(nx, dt, dx, UL, UR, UnL, UnR, D)