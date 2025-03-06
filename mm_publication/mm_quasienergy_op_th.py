import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
import time
#import seaborn as sns

start_time = time.time()

# Constants and parameters
w0 = 30000

# Domain for variable d
#z_list = [3.5, 10]
z_list = [0, 100]

num_list = range(0, 201)

Nx = 25
num_blocks = 20

# Phononic operators
ax = qt.tensor(qt.qeye(2), qt.destroy(Nx))

# Electronic states
sx = qt.tensor(qt.sigmax(), qt.qeye(Nx))
sy = qt.tensor(qt.sigmay(), qt.qeye(Nx))
sz = qt.tensor(qt.sigmaz(), qt.qeye(Nx))

# Function to generate the quasienergy operator matrix
def quasienergy_op_symmetric(num_blocks, block_matrices, eta):
    Qop0 = 0
    Qop1 = 0
    Qop2 = 0
    dim = 2 * num_blocks + 1

    for i in range(dim):
        Qop0 += qt.tensor((block_matrices[0] + (i - np.rint(0.5 * dim + 0.1))*eta), qt.basis(dim, i) * qt.basis(dim, i).dag())  # Diagonal terms

        next_index = (i + 1) % dim  # Periodic boundary condition
        Qop1 += qt.tensor(block_matrices[1], qt.basis(dim, next_index) * qt.basis(dim, i).dag())

        next_index_2 = (i + 2) % dim  # Periodic boundary for second-nearest neighbor
        Qop2 += qt.tensor(block_matrices[3], qt.basis(dim, next_index_2) * qt.basis(dim, i).dag())

    return Qop0 + Qop1 + Qop1.dag() + Qop2 + Qop2.dag()


total_iterations = len(z_list) * len(num_list)
iteration_count = 0

for g in z_list:
    for x in num_list:
        loop_start_time = time.time()
        
        eta = 3.5 + x / 25
        d = g
        # eta = g
        # d = x/2

        # Define Hamiltonians
        Hin0 = w0 * sz - (9 / 32) * d * sx
        Hex0 = ax.dag() * ax
        Hco0 = (1 / 4) * (ax.dag() + ax) * sx
        Hin1 = ((3 * eta * d) / 16) * sx
        Hex1 = -(3 / 8) * (ax.dag() + ax) * (ax.dag() - ax) + ((3 * d) / 64) * (ax.dag() - ax)
        Hco1 = -(3 / 2) * eta * (ax.dag() + ax) * sx
        Hex2 = -(9 / 64) * (ax.dag() + ax) * (ax.dag() + ax) + ((9 * d) / 256) * (ax.dag() + ax)

        block_matrices = [Hin0 + Hex0 + Hco0, (Hin1 + Hex1 + Hco1).dag(), Hin1 + Hex1 + Hco1, Hex2]
        Qop = quasienergy_op_symmetric(num_blocks, block_matrices, eta)

        # Compute eigenvalues and eigenstates
        eigvals, eigvecs = Qop.eigenstates()

        for n in range(Nx):
            phonon_state = qt.basis(Nx, n)
            nS_state = qt.tensor(qt.basis(2, 1), phonon_state, qt.basis(2 * num_blocks + 1, num_blocks + 1))
            nS_projector = nS_state*nS_state.dag()

            expectation_values = np.array([qt.expect(nS_projector, eigvecs[k]) for k in range(len(eigvals))])
            shifted_eigenvalues = eigvals + n
            # Save ordered overlap data
            overlap_filename = f'C:/Users/wsant/OneDrive/Dokumente/2024mmotion/1d_model/finite_temp/data_finite/overlap_d={d}_eta={eta}_n={n}'
            np.savetxt(overlap_filename, np.c_[shifted_eigenvalues, expectation_values])

        iteration_count += 1
        loop_end_time = time.time()
        print(f"Iteration {iteration_count}/{total_iterations} completed in {loop_end_time - loop_start_time:.2f} seconds")

total_time = time.time() - start_time
print(f"Total execution time: {total_time:.2f} seconds")





