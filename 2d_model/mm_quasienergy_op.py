import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
import time

start_time = time.time()

dx_list = range(0, 351)

#parameters of the trap: strength \gamma_i, defined as \omega_i/\omega 
gx, gy = 2, 2

#parameter \epsilon characterizing the asymmetry on the gradient
ep = 0

#harmonic oscillation length L = L_z, and characteristic electronic length \ell
hol, ell = 10, 100

#energy of states |nS_1> and |nP_1>
EnS1, EnP1 = 0, 1e5

#micromotion frequency \nu and trap frequency \omega
v, w = 10, 1

#number of bosons
Nx = 20
Ny = 20

#phononic operators
ax = qt.tensor(qt.qeye(2), qt.destroy(Nx), qt.qeye(Ny))
ay = qt.tensor(qt.qeye(2), qt.qeye(Nx), qt.destroy(Ny)) 

#electronic states are respectively defined through the set below 
sx = qt.tensor(qt.sigmax(), qt.qeye(Nx), qt.qeye(Ny)) 
sy = qt.tensor(qt.sigmay(), qt.qeye(Nx), qt.qeye(Ny))
sz = qt.tensor(qt.sigmaz(), qt.qeye(Nx), qt.qeye(Ny))

# Parameters
block_size = 2*Nx*Ny
num_blocks = 49  

def quasienergy_op(block_matrices):
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

for dx in dx_list:
    dy = 0
    #shorthand notations
    kx = dx/hol
    ky = dy/hol
    rho = ell/(3*hol)
    E = (EnS1 - EnP1)/2

    zx = np.sqrt(gx**2 + gy**2 + 1)/(2*gx)
    zy = np.sqrt(gx**2 + gy**2 + 1)/(2*gy)
    ux = (1 + ep)/(4*gx)
    uy = (1 - ep)/(4*gx)
    xi = v/w
    
    ############### internal motion #################
    Hin = E*sz + 0.25*rho*(kx*sx + ky*sy)

    ############### external motion #################
    Hex0 = gx*ax.dag()*ax + gy*ay.dag()*ay + (kx*ux*np.sqrt(gx))*(ax.dag() + ax) + (ky*uy*np.sqrt(gy))*(ay.dag() + ay)

    #time dependent on sine and cosine
    Hex1t = -1j*((zx*gx)*(ax.dag()*ax.dag() - ax*ax) - (zy*gy)*(ay.dag()*ay.dag() - ay*ay))
    Hex2t = 0.5*(zx**2)*gx*(ax.dag()*ax.dag() + 2*ax.dag()*ax + ax*ax) + 0.5*(zy**2)*gy*(ay.dag()*ay.dag() + 2*ay.dag()*ay + ay*ay)

    ############### coupling dynamics ###############
    Hco0 = -rho*(ux*np.sqrt(gx)*(ax.dag() + ax)*sx + uy*np.sqrt(gy)*(ay.dag() + ay)*sy) 

    #time dependent
    Hcot = rho*xi*((zx*np.sqrt(gx))*(ax.dag() + ax)*sx - (zy*np.sqrt(gy))*(ay.dag() + ay)*sy)
    
    #generating time dependent Hamiltonians
    H0 = Hin + Hex0 + Hco0

    # Define the block matrices
    M0 = H0        # Convert to Qobj
    Mp1 = Hcot - 1j*Hex1t
    Mm1 = Hcot + 1j*Hex1t
    M2 = Hex2t

    block_matrices = [M0, 0.5*Mp1, 0.5*Mm1, 0.5*M2]

    # Create the general matrix
    Qop = quasienergy_op(block_matrices)

    eigvl = Qop.eigenenergies()

    #state |nS> transition 
    nS_state = qt.basis(2*Nx*Ny*num_blocks, 0) 
    #projector unto |nS>    
    PnS = nS_state * nS_state.dag()

    # Compute expectation values and save to file
    ns_trans = []
    for k in range(len(eigvl)):
        ns_trans.append(qt.expect(PnS, Qop.eigenstates()[1][k]))

    #saving data
    nS_transition_txt = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\data\\nS_transition_dx={dx}_N={Nx}_Nb={num_blocks}_ep={ep}'
    np.savetxt(nS_transition_txt, np.c_[eigvl, ns_trans])
    print(f"--- {dx} done in {time.time() - start_time} seconds ---")

print("All calculations completed.")

