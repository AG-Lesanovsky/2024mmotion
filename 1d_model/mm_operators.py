#script with the set of operators for the main mm_dynamics
import qutip as qt
import numpy as np
from mm_parameters import*

#phononic operators
ax = qt.tensor(qt.qeye(2), qt.destroy(Nx))
 
#electronic states are respectively defined through the set bellow 
sx = qt.tensor(qt.sigmax(), qt.qeye(Nx)) 
sy = qt.tensor(qt.sigmay(), qt.qeye(Nx))
sz = qt.tensor(qt.sigmaz(), qt.qeye(Nx))

#Hamiltonians for internal, external motion and coupling
#Coupling and external Hamiltonian have time dependent and stationary part
#The former is denoted with Hxxt and the latter with Hxx0

############### internal motion #################
Hin = E*sz + 0.25*rho*kx*sx 

############### external motion #################
Hex0 = gx*ax.dag()*ax + kx*ux*np.sqrt(gx)*(ax.dag() + ax) 

#time dependent on sine and cosine
Hex1t = -1j*(zx*gx)*(ax.dag()*ax.dag() - ax*ax) 
Hex2t = 0.5*(zx**2)*gx*(ax.dag()*ax.dag() + 2*ax.dag()*ax + ax*ax) 

############### coupling dynamics ###############
Hco0 = -rho*ux*np.sqrt(gx)*(ax.dag() + ax)*sx 

#time dependent
Hcot = rho*xi*(zx*np.sqrt(gx))*(ax.dag() + ax)*sx 

#functions sine and cosine utilized in the code
sinvt, cosvt, cos2vt = 'np.sin({}*t)', 'np.cos({}*t)', 'np.cos({}*2*t)' 

sinvt, cosvt, cos2vt = sinvt.format(v), cosvt.format(v), cos2vt.format(v) 

#generating time dependent Hamiltonians
H0 = Hin + Hex0 + Hco0

H = [H0, [Hex1t, sinvt], [Hex2t, cos2vt], [Hcot, cosvt]]