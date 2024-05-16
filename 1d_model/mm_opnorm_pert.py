#Trace norm between the quasienergy operator and the unperturbed quasienergy operator

import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
import time

start_time = time.time()

dx_list = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000]
op_norm_values = []

#parameters of the trap: strength \gamma_i, defined as \omega_i/\omega 
gx = 2

#parameter \epsilon characterizing the assymmetry on the gradient
ep = 0

#harmonic oscillation length L = L_z, and characteristic electronic length \ell
hol, ell = 10, 100

#parameters \delta_x, \delta_y 
#dx= 100

#energy of states |nS_1> and |nP_1>
EnS1, EnP1 = 0, 1e5

#micromotion frequency \nu and trap frequency \omega
v, w = 10, 1

#number of bosons
Nx = 50

#phononic operators
ax = qt.tensor(qt.qeye(2), qt.destroy(Nx))
 
#electronic states are respectively defined through the set bellow 
sx = qt.tensor(qt.sigmax(), qt.qeye(Nx)) 
sy = qt.tensor(qt.sigmay(), qt.qeye(Nx))
sz = qt.tensor(qt.sigmaz(), qt.qeye(Nx))

for dx in dx_list:
    def quasienergy_op(block_size, num_blocks, block_matrices):
        matrix_list = []

        for i in range(num_blocks):
            row_list = []
            for j in range(num_blocks):
                if i == j:
                    block = block_matrices[0] + (np.rint(0.5*num_blocks + 0.1) - i)* (v/w)  # central block always has the identity times the corresponding frequency
                elif (i - j) == 1:
                    block = block_matrices[1]
                elif (i - j) == -1:
                    block = block_matrices[2]
                elif abs(i - j) == 2:
                    block = block_matrices[3]
                else:
                    block = np.zeros((block_size, block_size))
                row_list.append(block)

            matrix_list.append(row_list)

        return qt.Qobj(np.block(matrix_list))
    
    #shorthand notations
    kx = dx/hol
    rho = ell/(3*hol)
    E = (EnS1 - EnP1)/2

    zx = np.sqrt(gx**2 + 1)/(2*gx)
    ux = (1 + ep)/(4*gx)
    xi = v/w
    
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
    
    #generating time dependent Hamiltonians
    H0 = Hin + Hex0 + Hco0

    # Define the block matrices
    M0 = H0        
    Mp1 = Hcot - 1j*Hex1t
    Mm1 = Hcot + 1j*Hex1t
    M2 = Hex2t

    block_matrices = [M0, 0.5*Mp1, 0.5*Mm1, 0.5*M2]
    block_matrices_0 = [E*sz, 0*M0, 0*M0, 0*M0]
    # Parameters
    block_size = M0.shape[0]
    num_blocks = 49  # Example number of blocks

    # Quasienergy operator
    Qop = quasienergy_op(block_size, num_blocks, block_matrices)

    # Unperturbed quasienergy operator
    Qop0 = quasienergy_op(block_size, num_blocks, block_matrices_0)
    	
    op_norm = ((Qop - Qop0).norm())/(Qop0).norm()
    
    op_norm_values.append(op_norm)

op_norm_dx = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\data\\op_norm_dxvar'
np.savetxt(op_norm_dx, np.c_[dx_list, op_norm_values])
    
plt.figure(figsize=(7,5))

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24

#plt.ylim(0, 1)
plt.xlabel(r'$\delta_{x}$ (nm)', fontsize=24)
plt.ylabel(r'$\left\|\mathbb{V}\right\|/\left\|\mathbb{H}_\mathrm{F}\right\|$', fontsize=24)

# Plot the vertical lines
plt.plot(dx_list, op_norm_values, color='#b2182b', linestyle='solid', linewidth=1.5, marker='^', markersize=10)

#plt.xscale('log')

plt.tight_layout()
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))