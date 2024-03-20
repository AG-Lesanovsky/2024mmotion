#script with the set of operators for the main mm_dynamics
import qutip as qt
import numpy as np
from mm_parameters import*

#phononic operators
ax = qt.tensor(qt.qeye(2), qt.destroy(Nx))

#states |nS_1>, |nP_1> are respectively defined through the basis bellow
S1 = qt.tensor(qt.basis(2,0), qt.qeye(Nx)) 
P1 = qt.tensor(qt.basis(2,1), qt.qeye(Nx))

#dual vectors are defined subsequentely 
S1d = S1.dag()
P1d = P1.dag()

#Hamiltonians for internal, external motion and coupling
#Coupling and external Hamiltonian have time dependent and stationary part
#The former is denoted with Hxxt and the latter with Hxx0

############### internal motion #################
Hin = EnS1*S1*S1d + EnP1*P1*P1d
Hin = Hin + ((1 + ep)/4)*((dx*ell)/(hol**2))*(P1*S1d + S1*P1d)


############### external motion #################
Hex0 = dx*ax.dag()*ax + ((1 + ep)/(4*np.sqrt(gx)))*(dx/hol)*(ax.dag() + ax)

#time dependent on sine and cosine
Hex1t = -1j*(np.sqrt(1 + gx**2 + gy**2)/2)*(ax.dag()*ax.dag() - ax*ax)
Hex2t = -(1 + gx**2 + gx**2)*(1/(8*gx))*(ax.dag()*ax.dag() + 2*ax.dag()*ax + ax*ax) 

############### coupling dynamics ###############
Hco0 = -((1 + ep)/(4*np.sqrt(gx)))*(ell/hol)*(ax.dag() + ax)*(P1*S1d + S1*P1d)

#time dependent
Hcot = np.sqrt((1 + gx**2 + gy**2)/(4*gx))*((ell*v)/(hol*w))*(ax.dag() + ax)*(P1*S1d + S1*P1d)

#functions sine and cosine utilized in the code
sinvt, cosvt, cos2vt = 'np.sin({}*t)', 'np.cos({}*t)', 'np.cos({}*2*t)' 

sinvt, cosvt, cos2vt = sinvt.format(v), cosvt.format(v), cos2vt.format(v) 

#generating time dependent Hamiltonians
H0 = Hin + Hex0 + Hco0

H = [H0, [Hex1t, sinvt], [Hex2t, cos2vt], [Hcot, cosvt]]