import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
import time
from mm_operators import *
from mm_parameters import *

start_time = time.time()

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

# Define the block matrices
M0 = H0        # Convert to Qobj
Mp1 = Hcot - 1j*Hex1t
Mm1 = Hcot + 1j*Hex1t
M2 = Hex2t

block_matrices = [M0, 0.5*Mp1, 0.5*Mm1, 0.5*M2]

# Parameters
block_size = M0.shape[0]
num_blocks = 3  # Example number of blocks

# Create the general matrix
Qop = quasienergy_op(block_size, num_blocks, block_matrices)

eigvl, eigst = Qop.eigenstates()

#state |nS> transition 
nS_state = qt.basis(2*Nx*Ny*num_blocks, 0) 
#projector unto |nS>    
PnS = nS_state*nS_state.dag()

print(Qop)
#print(eigvl)
# Compute expectation values using vectorized operations
#ns_trans = np.abs(np.array([qt.expect(PnS, eigst[k]) for k in range(len(eigvl))])) ** 2
ns_trans = np.array([qt.expect(PnS, eigst[k]) for k in range(len(eigvl))])

################## Plotting of the transitions to the state nS #######################
plt.figure(figsize=(7,5))

plt.rcParams['mathtext.fontset'] = 'stix'
# Set the font family and size
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24

plt.ylim(0, 1)
#plt.xlim(-50000, -49900)

plt.xlabel(r'$\epsilon_{n} (\omega)$', fontsize=24)
plt.ylabel(r'$|\langle nS|u_{nm} \rangle|^2$', fontsize=24)

#print(sum(ns_trans))
print("--- %s seconds ---" % (time.time() - start_time))

plt.vlines(eigvl, 0, ns_trans, color='#b2182b', linestyle='solid', linewidth = 1.5)

plt.tight_layout()


plt.show()