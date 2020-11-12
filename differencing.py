#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 15:42:42 2020
@author: Nicholas Vieira
@differencing.py

Equations taken from PHYS643: Writing Hydro Codes (WHC) document. 
"""

import numpy as np


def __tridiag(ngrid, beta):
    return (1 + 2*beta)*np.eye(ngrid) - beta*(np.eye(ngrid, k=-1) + 
                                              np.eye(ngrid, k=1))    

def advec_FTCS(ngrid, dx, dt, u, f):
    # advance one time step using Forward-Time Central-Space method
    # equation (8)
    f[1:ngrid-1] = f[1:ngrid-1] - (0.5*u*dt/dx)*(f[2:] - f[:ngrid-2])
    return f


def advec_LF(ngrid, dx, dt, u, f):
    # propagate one time step using Lax-Friedrichs method
    # equation (11)
    f[1:ngrid-1] = 0.5*(f[2:] + f[:ngrid-2]) - (0.5*u*dt/dx)*(
            f[2:] - f[:ngrid-2])
    return f


def diffus_matr_build(ngrid, dx, dt, D):
    # compute beta term
    beta = D*dt/(dx**2)    
    # construct matrix A as described in section 2.4 of WHC document
    A = __tridiag(ngrid, beta)
    return A
    

def diffus(ngrid, A_inv, f):
    # solve for f^(n+1) given f^n according to the diffusion equation
    # equation (19)
    # the matrix is pre-built using diffus_matr_build() and has also been 
    # inverted
    f[1:ngrid-1] = np.matmul(A_inv[1:ngrid-1, 1:ngrid-1], f[1:ngrid-1])
    return f

