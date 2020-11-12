#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 16:53:17 2020
@author: Nicholas Vieira
@diffusion.py
"""

## imports
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc
rc('text', usetex=True)

## my script
from differencing import diffus_matr_build, diffus 

## CONSTANTS ##################################################################
NGRID = 50 # gridsize
NTSTEPS = 2500 # no. of time steps
DT = 1 # size of timestep
DX = 1 # size of spatial step
D_1 = 0.5 # diffusion coefficient #1
D_2 = 0.1 # diffusion coefficient #2
U = -0.1 # velocity


## FINITE DIFFERENCING SETUP ##################################################
# set up the 1D spatial grid
xgrid = np.arange(NGRID)*DX

# set the initial conditions
f1 = np.arange(NGRID)*DX/NGRID # for D=D_1
f2 = np.arange(NGRID)*DX/NGRID # for D=D_2

# build the matrices
A1 = diffus_matr_build(NGRID, DX, DT, D_1) # for D=D_1
A1[0][0] = 1 # set boundary conditions
A1[0][1] = 0
A1[-1][-1] = 1 + 2*D_1*DT/(DX**2)

A2 = diffus_matr_build(NGRID, DX, DT, D_2) # for D=D_1
A2[0][0] = 1 # set boundary conditions
A2[0][1] = 0
A2[-1][-1] = 1 + D_2*DT/(DX**2)

# invert the matrices
A1_inv = np.linalg.inv(A1)
A2_inv = np.linalg.inv(A2)


## PLOTTING ###################################################################
# set up the plots
plt.ion() # interactive on
fig, axes = plt.subplots(1,2, figsize=(12,12)) # 1 x 2 array of subplots

# initial conditions as faded grey lines, for reference
axes[0].plot(xgrid, f1, color="k", ls="-", alpha=0.7)
axes[1].plot(xgrid, f2, color="k", ls="-", alpha=0.7)

# fix y limits for easy comparison
axes[0].set_ylim(-0.05, 1.05)
axes[1].set_ylim(-0.05, 1.05)

# plotting objects which will be updated at each timestep
plt1, = axes[0].plot(xgrid, f1, color="#0165fc", marker="o")
plt2, = axes[1].plot(xgrid, f2, color="#0165fc", marker="o")

# allow animations
fig.canvas.draw()

# give the subplots titles
axes[0].set_title(r"$D = "+f"{D_1}"+r"$", fontsize=20)
axes[1].set_title(r"$D = "+f"{D_2}"+r"$", fontsize=20)

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
        ticker.MultipleLocator(NGRID//5)) # enforce x-axis has exactly 6 ticks
axes[1].xaxis.set_major_locator(
        ticker.MultipleLocator(NGRID//5))


## TIME EVOLUTION #############################################################
for i in range(NTSTEPS):
    
    # update f1 (D=D_1) with implicit method for diffusion
    f1 = diffus(NGRID, A1_inv, f1)
    #f1 = diffus(NGRID, A1, f1)
    
    # update f2 (D=D_2) with implicit method for diffusion
    f2 = diffus(NGRID, A2_inv, f2)
    #f2 = diffus(NGRID, A2, f2)
    
    # update the plots
    plt1.set_ydata(f1)
    plt2.set_ydata(f2)
    
    fig.canvas.draw() # draw it
    plt.pause(0.005)