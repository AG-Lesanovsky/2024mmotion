#script with the set of parameters for the main mm_dynamics
import qutip as qt
import numpy as np

#parameters of the trap: strength \gamma_i, defined as \omega_i/\omega 
gx = 2

#parameter \epsilon characterizing the assymmetry on the gradient
ep = 0

#harmonic oscillation length L = L_z, and characteristic electronic length \ell
hol, ell = 10, 100

#parameters \delta_x, \delta_y 
dx= 100

#energy of states |nS_1> and |nP_1>
EnS1, EnP1 = 0, 1e5

#micromotion frequency \nu and trap frequency \omega
v, w = 10, 1

#number of bosons
Nx = 50

#shorthand notations
kx = dx/hol
rho = ell/(3*hol)
E = (EnS1 - EnP1)/2

zx = np.sqrt(gx**2 + 1)/(2*gx)
ux = (1 + ep)/(4*gx)
xi = v/w