#script with the set of parameters for the main mm_dynamics
import qutip as qt
import numpy as np
from mm_parameters import*

#phononic operators
ax = qt.tensor(qt.qeye(4), qt.destroy(Nx), qt.qeye(Ny), qt.qeye(Nz))
ay = qt.tensor(qt.qeye(4), qt.qeye(Nx), qt.destroy(Ny), qt.qeye(Nz)) 
az = qt.tensor(qt.qeye(4), qt.qeye(Nx), qt.qeye(Ny), qt.destroy(Nz))

#states |nS_1, -1>, |nS_1, 1>, |nP_1, -1> and |nP_1, 1> are respectively defined through the basis bellow
S1m = qt.tensor(qt.basis(4,0), qt.qeye(Nx), qt.qeye(Ny), qt.qeye(Nz)) 
S1p = qt.tensor(qt.basis(4,1), qt.qeye(Nx), qt.qeye(Ny), qt.qeye(Nz))
P1m = qt.tensor(qt.basis(4,2), qt.qeye(Nx), qt.qeye(Ny), qt.qeye(Nz))
P1p = qt.tensor(qt.basis(4,3), qt.qeye(Nx), qt.qeye(Ny), qt.qeye(Nz))

#dual vectors are defined subsequentely 
S1md = S1m.dag()
S1pd = S1p.dag()
P1md = P1m.dag()
P1pd = P1p.dag()

#Hamiltonians for internal, external motion and coupling
#Coupling and external Hamiltonian have time dependent and stationary part
#The former is denoted with Hxxt and the latter with Hxx0

############### internal motion #################
Hin = EnS1*(S1m*S1md + S1p*S1pd) + EnP1*(P1m*P1md + P1p*P1pd)
Hin = Hin + ((1 + ep)/4)*((dx*ell)/(hol**2))*(P1p*S1md + P1m*S1pd + S1p*P1md + S1m*P1pd)
Hin = Hin + ((1 + ep)/4)*((dx*ell)/(hol**2))*(P1p*S1md - P1m*S1pd + S1p*P1md + S1m*P1pd)


############### external motion #################
Hex0 = dx*ax.dag()*ax + dy*ay.dag()*ay + az.dag()*az + ((1 + ep)/(4*np.sqrt(gx)))*(dx/hol)*(ax.dag() + ax) + ((1 - ep)/(4*np.sqrt(gy)))*(dy/hol)*(ay.dag() + ay)

#time dependent on sine and cosine
Hex1t = 1j*(np.sqrt(1 + gx**2 + gy**2)/2)*(ay.dag()*ay.dag() - ay*ay- ax.dag()*ax.dag() + ax*ax)
Hex2t = -(1 + gx**2 + gy**2)*(1/(8*gx)*(ax.dag()*ax.dag() + 2*ax.dag()*ax + ax*ax) + 1/(8*gy)*(ay.dag()*ay.dag() + 2*ay.dag()*ay + ay*ay))

############### coupling dynamics ###############
Hco0 = -((1 + ep)/(4*np.sqrt(gx)))*(ell/hol)*(ax.dag() + ax)*(P1p*S1md + P1m*S1pd + S1p*P1md + S1m*P1pd)
Hco0 = Hco0 - ((1 - ep)/(4*np.sqrt(gy)))*(ell/hol)*(ay.dag() + ay)*(P1p*S1md - P1m*S1pd + S1p*P1md - S1m*P1pd)
Hco0 = Hco0 - 0.5*(ell/hol)*(az.dag() + az)*(P1m*S1md - P1p*S1pd + S1m*P1md - S1p*P1pd)

#time dependent
Hcot = (ax.dag() + ax)*(P1p*S1md + P1m*S1pd + S1p*P1md + S1m*P1pd)
Hcot = Hcot - 1j*(ay.dag() + ay)*(P1p*S1md - P1m*S1pd + S1p*P1md - S1m*P1pd)
Hcot = np.sqrt((1 + gx**2 + gy**2)/(4*gx))*((ell*v)/(hol*w))*Hcot

#functions sine and cosine utilized in the code
sinvt, cosvt, cos2vt = 'np.sin({}*t)', 'np.cos({}*t)', 'np.cos({}*2*t)' 

sinvt, cosvt, cos2vt = sinvt.format(v), cosvt.format(v), cos2vt.format(v) 

#generating time dependent Hamiltonians
H0 = Hin + Hex0 + Hco0

H = [H0, [Hex1t, sinvt], [Hex2t, cos2vt], [Hcot, cosvt]]