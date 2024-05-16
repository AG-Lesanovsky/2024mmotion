import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
import time

start_time = time.time()

#dx_list = [301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350]
dx_list = range(0, 351)

#parameters of the trap: strength \gamma_i, defined as \omega_i/\omega 
gx = 2

#parameter \epsilon characterizing the assymmetry on the gradient
ep = 0.5

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
    M0 = H0        # Convert to Qobj
    Mp1 = Hcot - 1j*Hex1t
    Mm1 = Hcot + 1j*Hex1t
    M2 = Hex2t

    block_matrices = [M0, 0.5*Mp1, 0.5*Mm1, 0.5*M2]

    # Parameters
    block_size = M0.shape[0]
    num_blocks = 49  # Example number of blocks

    # Create the general matrix
    Qop = quasienergy_op(block_size, num_blocks, block_matrices)

    eigvl, eigst = Qop.eigenstates()

    #state |nS> transition 
    nS_state = qt.basis(2*Nx*num_blocks, 0) 
    #projector unto |nS>    
    PnS = nS_state*nS_state.dag()

    print(Qop)
    #print(eigvl)
    # Compute expectation values using vectorized operations
    #ns_trans = np.abs(np.array([qt.expect(PnS, eigst[k]) for k in range(len(eigvl))])) ** 2
    ns_trans = np.array([qt.expect(PnS, eigst[k]) for k in range(len(eigvl))])

    #print(sum(ns_trans))
    print("--- %s seconds ---" % (time.time() - start_time))

    #saving data
    nS_transition_txt = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\data\\nS_transition_dx={dx}_N={Nx}_Nb={num_blocks}_ep={ep}'
    np.savetxt(nS_transition_txt, np.c_[eigvl, ns_trans])


#plt.show()