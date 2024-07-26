import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
import time

start_time = time.time()

#This is the most important code. It generates all data used in our problem, diagonalizing the matrix associated with the quasienergy operator.
#Each calculation is optmized and detailed.


num_list = range(0, 201)

w0 = 30000
#z_list = [0, 50, 100]
z_list = [10]

# Number of bosons
Nx = 50
num_blocks = 35  # Half the number of blocks for symmetric range

# Phononic operators
ax = qt.tensor(qt.qeye(2), qt.destroy(Nx))

# Electronic states are respectively defined through the set below 
sx = qt.tensor(qt.sigmax(), qt.qeye(Nx)) 
sy = qt.tensor(qt.sigmay(), qt.qeye(Nx))
sz = qt.tensor(qt.sigmaz(), qt.qeye(Nx))

def quasienergy_op_symmetric(block_size, num_blocks, block_matrices, eta):
    total_blocks = 2 * num_blocks + 1
    matrix_list = []

    for i in range(-num_blocks, num_blocks + 1):
        row_list = []
        for j in range(-num_blocks, num_blocks + 1):
            if i == j:
                block = block_matrices[0] + (np.rint(0.5 * total_blocks + 0.1) - i) * eta  # central block always has the identity times the corresponding frequency
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

for g in z_list:
    for x in num_list:
        #eta = 3.5 + x / 25
        eta = g
        #d = g
        d = x/2

        ############### Zero mode #################
        Hin0 = w0 * sz - (9 / 32) * d * sx
        Hex0 = ax.dag() * ax
        Hco0 = (1 / 4) * (ax.dag() + ax) * sx

        ############### First mode #################
        Hin1 = ((3*eta*d)/16)*sx
        Hex1 = -(3 / 8) * (ax.dag() + ax) * (ax.dag() - ax) +  ((3*d)/64) * (ax.dag() - ax) 
        Hco1 = -(3 / 2) * eta * (ax.dag() + ax) * sx

        ############### Second mode ###############
        Hex2 = -(9 / 64) * (ax.dag() + ax) * (ax.dag() + ax) + ((9*d)/256)*(ax.dag() + ax)

        # Define the block matrices
        M0 = Hin0 + Hex0 + Hco0
        Mm1 = Hin1 + Hex1 + Hco1
        Mp1 = Mm1.dag()
        M2 = Hex2

        block_matrices = [M0, Mp1, Mm1, M2]

        # Parameters
        block_size = M0.shape[0]

        # Create the general matrix
        Qop = quasienergy_op_symmetric(block_size, num_blocks, block_matrices, eta)

        eigvals, eigvecs = Qop.eigenstates()

        # State |nS> transition 
        nS_state = qt.tensor(qt.basis(2, 1), qt.basis(Nx), qt.basis(2 * num_blocks + 1))

        # Set dimensions for projector
        nS_state.dims = [[2 * Nx * (2 * num_blocks + 1)], [1]]

        # Projector onto |nS>
        PnS = nS_state * nS_state.dag()

        # Compute expectation values using vectorized operations
        ns_trans = np.array([qt.expect(PnS, eigvecs[k]) for k in range(len(eigvals))])

        print(Qop)
        print(Qop.tr())
        print("--- %s seconds ---" % (time.time() - start_time))

        # Saving data
        nS_transition_txt = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\data\\smnS_transition_d={d}_eta={eta}'
        np.savetxt(nS_transition_txt, np.c_[eigvals, ns_trans])


