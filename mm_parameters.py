#script with the set of parameters for the main mm_dynamics
import qutip as qt
import numpy as np

#parameters of the trap: strength \gamma_i, defined as \omega_i/\omega 
gx, gy, gz = 2, 2, 1

#parameter \epsilon characterizing the assymmetry on the gradient
ep = 0

#harmonic oscillation length L = L_z, and characteristic electronic length \ell
hol, ell = 10, -20

#parameters \delta_x, \delta_y 
dx, dy = 0, 0

#energy of states |nS_1> and |nP_1>
EnS1, EnP1 = 0, 1e5

#micromotion frequency \nu and trap frequency \omega
v, w = 10, 1

#number of bosons for x, y, and z direction
Nx, Ny, Nz = 10, 10, 10

#shorthand notations
kx, ky = dx/hol, dy/hol
rho = ell/(3*hol)
E = (EnS1 - EnP1)/2

zx, zy, zz = np.sqrt(gx**2 + gy**2 + 1)/(2*gx), np.sqrt(gx**2 + gy**2 + 1)/(2*gy), np.sqrt(gx**2 + gy**2 + 1)/2 
ux, uy, uz = (1 + ep)/(4*gx), (1 - ep)/(4*gx), 0.25
xi = v/w