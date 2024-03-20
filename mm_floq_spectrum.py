#script to run Floquet eigenspectrum of the time-dependent Hamiltonian   
#everything is in units of hbar times trap frequency \omega
import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
import time
#from mm_operators import*
from mm_operators1d import*
from mm_parameters import*

start_time = time.time()

#calculating the Floquet quasi-energies and Floquet modes associated to the system
dx_vec = np.linspace(0, 1, 50)

hs_dim = Nx*2

q_energies = np.zeros((len(dx_vec), hs_dim))
T = 2*np.pi/v 

for idx, dx in enumerate(dx_vec):

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

    opts = qt.Options(nsteps = 500000)
    f_modes, f_energies = qt.floquet_modes(H, T, options=opts)  
    q_energies[idx,:] = f_energies


plt.figure() 
for k in range(0, hs_dim):
    plt.plot(dx_vec, (((q_energies[:,k]/v) + 0.15) % 1) - 0.15, color='black', marker='.', linestyle='') 
plt.xlabel(r'$\delta_x$') 
plt.ylabel(r'$E$') 
plt.title(r'Floquet quasienergies')

print("--- %s seconds ---" % (time.time() - start_time))

plt.show() 