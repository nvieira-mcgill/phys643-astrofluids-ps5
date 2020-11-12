#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 15:42:42 2020
@author: Nicholas Vieira
@differencing.py

Functions which are called by advection.py and diffusion.py. Functions perform
various steps of finite differencing, including:
    - Forward-Time Central-Space method, for advection
    - Lax-Friedrichs method, for advection
    - An implicit method for solving the diffusion equation

All equations are taken from the "PHYS643: Writing Hydro Codes" (WHC) document.
The sections alluded to are those in the WHC document as well. 
"""

## imports
import numpy as np


def __tridiag(ngrid, beta):
    """
    ngrid: dimension of square matrix (int)
    beta:  D*dt/(dx)**2 , a parameter of relevance when solving the diffusion
           equation (float)
    
    Utility function which returns a tridiagonal matrix with value (-beta) 
    along the main diagonal and (1+2*beta) along the diagonals just above and 
    below the main diagonal.
    
    Matrix is constructed using the code snippet provided in section 2.5.2. 
    """
    
    return (1 + 2*beta)*np.eye(ngrid) - beta*(np.eye(ngrid, k=-1) + 
                                              np.eye(ngrid, k=1))
    

def advec_FTCS(ngrid, dx, dt, u, f):
    """
    ngrid: dimension of square matrix (int)
    dx:    spatial step size (float)
    dt:    temporal step size (float)
    u:     velocity (float)
    f:     the array to update
    
    Advance one time step dt using the Forward-Time Central-Space method for 
    finite differencing to solve the advection equation (Eqn. (8)). The 
    quantity f is changed in-place and then returned. 
    """    
    # advance one time step using Forward-Time Central-Space method
    f[1:ngrid-1] = f[1:ngrid-1] - (0.5*u*dt/dx)*(f[2:] - f[:ngrid-2])
    return f


def advec_LF(ngrid, dx, dt, u, f):
    """
    ngrid: dimension of square matrix (int)
    dx:    spatial step size (float)
    dt:    temporal step size (float)
    u:     velocity (float)
    f:     the array to update (array)
    
    Advance one time step dt using the Lax-Friedrichs method for 
    finite differencing to solve the advection equation (Eqn. (11)). The 
    quantity f is changed in-place and then returned. 
    """
    # propagate one time step using Lax-Friedrichs method
    f[1:ngrid-1] = 0.5*(f[2:] + f[:ngrid-2]) - (0.5*u*dt/dx)*(
            f[2:] - f[:ngrid-2])
    return f


def diffus_matr_build(ngrid, dx, dt, D):
    """
    ngrid: dimension of square matrix (int)
    dx:    spatial step size (float)
    dt:    temporal step size (float)
    D:     diffusion coefficient (float)
    
    Construct the matrix A which satisfies f**(n+1) = inverse(A) * f**n. Calls
    the utility function __tridiag() to do this. Used to solve the diffusion
    equation using an implicit method (Eqn. (19)). 
    """
    # compute beta term
    beta = D*dt/(dx**2)    
    # construct matrix A as described in section 2.4 of WHC document
    A = __tridiag(ngrid, beta)
    return A
    

def diffus(ngrid, A_inv, f):
    """
    ngrid: dimension of square matrix (int)
    A_inv: inverted A matrix, where A is the matrix which satisfies f**(n+1) = 
           A_inv * f**n in the implicit method for solving the diffusion 
           equation 
    f:     the array to update (array)
    
    Compute f**(n+1) given all the parameters of the space/time-grid encoded in
    the matrix A_inv which satisfies f**(n+1) = A_inv * f**n (Eqn. (19)). Used
    to solve the diffusion equation. The quantity f is changed in-place and 
    then returned.  
    """
    # solve for f^(n+1) given f^n according to the diffusion equation
    # matrix has already been built + inverted outside this function
    #f[1:ngrid-1] = np.matmul(A_inv[1:ngrid-1, 1:ngrid-1], f[1:ngrid-1])
    f[1:ngrid-1] = np.matmul(A_inv, f)[1:ngrid-1]
    return f


## ALTERNATIVE METHOD FOR SOLVING FOR f WHICH IS PROBABLY SLOWER BECAUSE YOU 
## KEEP RE-INVERTING MATRIX A, WHEREAS FUNCTION ABOVE USES PRE-INVERTED MATRIX
#def diffus(ngrid, A, f):
#    """
#    ngrid: dimension of square matrix (int)
#    A:     matrix A which satisfies f**(n+1) = inverse(A) * f**n in the 
#           implicit method for solving the diffusion equation 
#    f:     the array to update (array)
#    
#    Compute f**(n+1) given all the parameters of the space/time-grid encoded in
#    the matrix A_inv which satisfies f**(n+1) = A_inv * f**n (Eqn. (19)). Used
#    to solve the diffusion equation. The quantity f is changed in-place and 
#    then returned.  
#    """
#    # solve for f^(n+1) given f^n according to the diffusion equation
#    # matrix has already been built + inverted outside this function
#    f = np.linalg.solve(A, f)
#    return f
