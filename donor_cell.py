#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 22:01:41 2020
@author: Nicholas Vieira
@donor_cell.py

Functions for carrying out the donor cell advection scheme, a type of finite
volume method, to create a simple 1-dimensional hydrodynamics code solver. 
Used by the script Q5_1Dhydro.py.

All equations are taken from the "PHYS643: Writing Hydro Codes" (WHC) document.
The sections alluded to are those in the WHC document as well. 
"""

import numpy as np    
    
def J_compute(ngrid, u, f):
    """
    ngrid: dimension of square matrix (int)
    u:     velocity (float)
    f:     quantity of interest (array)
    
    Compute the flux J into a cell, as given in Section 3.1.1.
    """
    # initial: no flux
    J = np.zeros(ngrid-1)
    
    # where velocity u is POSITIVE (rightwards)
    idx_pos = np.where(u > 0)[0]
    J[idx_pos] = u[idx_pos]*f[idx_pos]
    
    # where velocity u is NEGATIVE (leftwards)
    idx_neg = np.where(u < 0)[0]
    J[idx_neg] = u[idx_neg]*f[idx_neg+1]
    
    return J  


def f_update_nosource(ngrid, dx, dt, u, f):
    """
    ngrid: dimension of square matrix (int)
    dx:    spatial step size (float)
    dt:    temporal step size (float)
    u:     velocity (float)
    f:     quantity of interest (array)
    
    Update quantity f in the donor-cell advection scheme, WITHOUT a source. f 
    may be f1 (rho) or f2 (rho * u). Following Equation (25).
    """
    # compute flux
    j = J_compute(ngrid, u, f)
    
    # update 
    f[1:ngrid-1] = f[1:ngrid-1] - (dt/dx)*(j[1:ngrid-1] - j[:ngrid-2])
    
    return f


def f_update_withsource(ngrid, dx, dt, u, f1, f2, cs2):
    """
    ngrid: dimension of square matrix (int)
    dx:    spatial step size (float)
    dt:    temporal step size (float)
    u:     velocity (float)
    f1:    first quantity of interest (array)
    f2:    second quantity of interest (array)
    cs2:   speed of sound in the medium, squared (float)
    
    Update quantity f2 in the donor-cell advection scheme, WITH a source. f2 
    is the quantity rho * u. Following Equation (28).
    
    Assumes that the two arrays have already been updated without a source 
    using the function f_update_nosource() above. 
    """     
    # update f2 using f1
    f2[1:ngrid-1] = f2[1:ngrid-1] - (dt*cs2/dx)*(f1[2:] - f1[:ngrid-2])
    return f1, f2
    

def gaussian_perturb(xgrid, mu, sigma, amp=1.0):
    """
    xgrid:  the spatial grid on which to define a Gaussian perturbation (array)
    mu:     mean of the Gaussian (float)
    sigma:  standard deviation of the Gaussian (float)
    amp:    constant amplitude applied to the Gaussian (float, default 1.0)
    
    Construct a Gaussian perturbation to some quantity defined on a spatial 
    grid given by xgrid.
    """
    # build + return the perturbation
    perturb = amp*np.exp(-((xgrid-mu)/sigma)**2)
    return perturb