#script with the set of parameters for the main mm_dynamics
import qutip as qt
import numpy as np

#parameters of the trap: strength \gamma_i, defined as \omega_i/\omega 
gx, gy, gz = 1, 1, 1

#parameter \epsilon characterizing the assymmetry on the gradient
ep = 0.1

#parameters \delta_x, \deltay 
dx, dy = 1, 1

#harmonic oscillation length L = L_z, and characteristic electronic length \ell
hol, ell = 1, 1

#energy of states |nS_1> and |nP_1>
EnS1, EnP1 = 1, 1

#micromotion frequency \nu and trap frequency \omega
v, w = 1, 1

#number of bosons for x, y, and z direction
Nx, Ny, Nz = 10, 10, 10
