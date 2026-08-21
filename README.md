FEM1D_padapting is a 1-dimensional FEM code for a very general wave equation problem. While the current form is emphasizing solutions to the scattering matrix equation problem, eventually it will be updated to 
solve the standing wave mode solutions to the eigenvalue matrix equation. The code is intended only as a learning experience for work in Notaros lab.

Features: 
1. Code can solve Ax = b  or Ax = lambda x.
2.  Most importantly the code will implement the hp adaptive refinement methods detailed by Cam Key, Jake Harmon, and Notaros is various papers.

It is ultimately a work in progress.

To run, simply set the length of your "wire" domain (L), the number of FEM elements you would like (M), and then the initial polynomial orders (n). Then choose either "scattering" or "waveguide" to solve which 
problem you desire. And then "HW" solution simply solves the scattering problem as described in ECE 540 assignment, but with adaptive anisotropic mesh refinement.
