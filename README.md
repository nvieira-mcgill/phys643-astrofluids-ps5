

**Nicholas Vieira**

**Python 3.7**

**Collaborator: Sabrina Berger**

================================================

This repository contains the following scripts which constitute a library of
functions:
- ```differencing.py```
- ```donor_cell.py```

These scripts are called by the following 3 scripts:
- ```Q3_advection.py```
- ```Q4_diffusion_advection.py```
- ```Q5_1Dhydro.py```

The first of which provides the solution to question #3 of the assignment, and so 
forth. 

All scripts can be run via the command line via e.g. 

	>>> python Q3_advection.py

With the appropriate dependencies, i.e., 
- matplotlib                3.0.3
- numpy                     1.16.2


**[5] 1-Dimensional Hydro solver**

As the amplitude of the Gaussian perturbation increases, the solution becomes 
increasingly unstable. In particular:
- For a density with baseline offset ```OFFSET=100.0``` and perturbation 
  amplitude ```AMP=0.5```, the solution is stable.
- For a density with the same base offset and a perturbation amplitude of 
  ```AMP=10.0```, the solution blows up, and the amplitude of the sound waves 
  tends to infinity. The hydro solver eventually crashes under these conditions.

A shock begins to appear depending on the ratio between dx and dt. In the scheme 
we used, no viscosity is explicitly included, but we *have* incurred a 
**"numerical viscosity"** due to the approximations used in the finite differencing
steps. This numerical viscosity is proportional to dx**2 / dt, and so it is this 
quantity which sets the width of the shock (Equation (15)). 
