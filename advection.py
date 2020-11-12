#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 12:08:17 2020
@author: Nicholas Vieira
@advection.py
"""

import numpy as np

from differencing import advec_FTCS, advec_LF

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from matplotlib import rc
rc('text', usetex=True)

## CONSTANTS
NGRID = 50 # gridsize
NTSTEPS = 2500 # no. of time steps
DT = 1 # size of timestep
DX = 1 # size of spatial step
U = -0.1 # velocity

## FINITE DIFFERENCING
# set up the 1D spatial grid
xgrid = np.arange(NGRID)*DX

# set the initial conditions
f1 = np.arange(NGRID)*DX/NGRID # for FTCS
f2 = np.arange(NGRID)*DX/NGRID # for Lax-Friedrich

# set up the plots
plt.ion() # interactive on
fig, axes = plt.subplots(1,2, figsize=(12,12)) # 1 x 2 array of subplots

# initial conditions as faded grey lines, for reference
axes[0].plot(xgrid, f1, color="k", ls="-", alpha=0.7)
axes[1].plot(xgrid, f2, color="k", ls="-", alpha=0.7)

# fix y limits for easy comparison
axes[0].set_ylim(-0.05, 2.05)
axes[1].set_ylim(-0.05, 2.05)

# plotting objects which will be updated at each timestep
plt1, = axes[0].plot(xgrid, f1, color="#ff474c", marker="o")
plt2, = axes[1].plot(xgrid, f2, color="#ff474c", marker="o")

# allow animations
fig.canvas.draw()

# give the subplots titles
axes[0].set_title("Forward-Time Central-Space", fontsize=20)
axes[1].set_title("Lax-Friedrichs", fontsize=20)

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

# time evolution
for i in range(NTSTEPS):
    
    # update f1 with FTCS
    f1 = advec_FTCS(NGRID, DX, DT, U, f1)
    
    # update f2 with Lax-Friedrich
    f2 = advec_LF(NGRID, DX, DT, U, f2)
    
    # update the plots
    plt1.set_ydata(f1)
    plt2.set_ydata(f2)
    
    fig.canvas.draw() # draw it
    plt.pause(0.005)




