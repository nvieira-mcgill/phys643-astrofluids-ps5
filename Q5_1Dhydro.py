#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 23:02:21 2020
@author: Nicholas Vieira
@Q5_1Dhydro.py

In this script, WHC = the "PHYS643: Writing Hydro Codes" document. All 
equations come from this document.

Using the donor cell advection scheme, carries out the propagation of a sound
wave in a uniform density medium. The sound wave is created by a small Gaussian
perturbation to the density. The quantities rho (density) and rho * u, where u
is the velocity, are animared.

This script provides an answer to Question 5 of the computing assignment for 
PHYS643.

All relevant constants (grid size, step size, etc.) are set at the beginning 
of this script.
"""

## imports
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc
rc('text', usetex=True)

## my script
from donor_cell import (J_compute,
                        f_update_nosource, f_update_withsource, 
                        gaussian_perturb)

## CONSTANTS ##################################################################
NGRID = 10000 # gridsize
NTSTEPS = int(1e5) # no. of time steps
DT = 1e-3 # size of timestep
DX = 10e-1 # size of spatial step
CS = 300 # speed of sound in medium
AMP = 1.0 # amplitude of the Gaussian perturbation
OFFSET = 100 # offset of the Gaussian perturbation

## SETUP ######################################################################
# set up the 1D spatial grid
xgrid = np.arange(NGRID)*DX

# set the initial conditions
# f1 is density rho
# induce a Gaussian perturbation on the density, where the density is given
# by some constant <OFFSET>
f1 = AMP*gaussian_perturb(xgrid, 
                          mu=np.mean(xgrid), 
                          sigma=np.std(xgrid)/10, 
                          amp=AMP) # the Gaussian perturbation
f1 += OFFSET # constant  

# set up initial velocity as all 0
u = np.zeros(NGRID-1)

# f2 is density rho * velocity u (all 0 as well, initially)
f2 = np.zeros(NGRID)
f2[:-1] = u*f1[:-1]


## PLOTTING ###################################################################
# set up the plots
plt.ion() # interactive on
fig, axes = plt.subplots(1,2, figsize=(12,12)) # 1 x 2 array of subplots

# initial conditions as faded grey lines, for reference
axes[0].plot(xgrid, f1, color="k", ls="-", alpha=0.5)
axes[1].plot(xgrid, f2, color="k", ls="-", alpha=0.5)

# plotting objects which will be updated at each timestep
plt1, = axes[0].plot(xgrid, f1, color="#04d8b2", marker="o", markersize=1, 
            lw=1)
plt2, = axes[1].plot(xgrid, f2, color="#04d8b2", marker="o", markersize=1,
            lw=1)

# allow animations
fig.canvas.draw()

# give the subplots titles
axes[0].set_title(r"$f_1 = \rho$", fontsize=20)
axes[1].set_title(r"$f_2 = \rho u$", fontsize=20)

# misc pretty plotting
axes[0].grid() # add a grid to both plots
axes[1].grid()
axes[0].set_xlabel(r"$x$", fontsize=16) # set x-axis label
axes[1].set_xlabel(r"$x$", fontsize=16)
axes[0].set_ylabel(r"$f$", fontsize=16) # set y-axis label (only for left plot)
axes[0].tick_params(axis="x", labelsize=20) # bigger tick labels
axes[0].tick_params(axis="y", labelsize=20)
axes[1].tick_params(axis="x", labelsize=20)
axes[1].tick_params(axis="y", labelsize=20)
axes[0].xaxis.set_major_locator(
        ticker.MultipleLocator((xgrid[-1]+1)/5)) 
axes[1].xaxis.set_major_locator(
        ticker.MultipleLocator((xgrid[-1]+1)/5))


## TIME EVOLUTION #############################################################
minf2 = min(f2) # initial value for minimum of f2
maxf2 = max(f2) # initial value for maixmum of f2
for i in range(NTSTEPS):

    # compute velocity array according to Eqn. (26)
    u[:NGRID-1] = 0.5 * ((f2[:NGRID-1]/f1[:NGRID-1]) + 
                         (f2[1:NGRID]/f1[1:NGRID]))    

    # compute fluxes
    j1 = J_compute(NGRID, u, f1)
    j2 = J_compute(NGRID, u, f2)

    # update f_1, f_2 WITHOUT source term
    # according to Equations (25), (26)
    f1 = f_update_nosource(ngrid=NGRID, dx=DX, dt=DT, u=u, f=f1)
    f2 = f_update_nosource(ngrid=NGRID, dx=DX, dt=DT, u=u, f=f2)
    
    # update f_1, f_2 WITH source term
    # according to Equation (28)
    f1, f2 = f_update_withsource(ngrid=NGRID, dx=DX, dt=DT, u=u, f1=f1, f2=f2,
                                 cs2=CS**2)
    
    # enforce reflective boundary conditions as described in Section 3.1.3 
    f1[0] = f1[0] - (DT/DX)*j1[0]
    f1[-1] = f1[-1] + (DT/DX)*j1[-1]
    f2[0] = f2[0] - (DT/DX)*j2[0]
    f2[-1] = f2[-1] + (DT/DX)*j2[-1]
    
    # update the plots
    if i % 50 == 0: # on every 50th iter
        plt1.set_ydata(f1)
        plt2.set_ydata(f2)
        # set the minimum and maximum of the y-axis for f2 to the minimum 
        # and maximum across ALL time steps
        minf2 = min(f2.tolist()+[minf2])
        maxf2 = max(f2.tolist()+[maxf2])
        #axes[1].set_ylim(1.1*minf2, 1.1*maxf2)
        fig.canvas.draw() # draw it
        plt.pause(0.001)