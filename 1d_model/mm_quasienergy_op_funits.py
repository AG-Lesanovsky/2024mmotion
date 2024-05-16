import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
import time

start_time = time.time()

d_list = range(0, 351)

w, mu, v, g = 1, 0.5, 5, 9/50
wsp = 25e3

#number of bosons
Nx = 50

#phononic operators
ax = qt.tensor(qt.qeye(2), qt.destroy(Nx))
 
#electronic states are respectively defined through the set bellow 
sx = qt.tensor(qt.sigmax(), qt.qeye(Nx)) 
sy = qt.tensor(qt.sigmay(), qt.qeye(Nx))
sz = qt.tensor(qt.sigmaz(), qt.qeye(Nx))

for d in d_list:
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
    
    ############### internal motion #################
    Hin = wsp*sz + 3*((mu**2)/g)*d*sx 

    ############### external motion #################
    Hex0 = w*ax.dag()*ax + ((3*d)/2)*np.sqrt((mu**2)/(g*w))*(ax.dag() + ax) 

    #time dependent on sine and cosine
    Hex1t = -0.5j*np.sqrt(w**2 + mu**2)*(ax.dag()*ax.dag() - ax*ax) 
    
    Hex2t = -0.5*((w**2 + mu**2)/(4*w))*(ax.dag()*ax.dag() + 2*ax.dag()*ax + ax*ax) 

    ############### coupling dynamics ###############
    Hco0 = -mu*np.sqrt((mu**2)/(4*g*w))*(ax.dag() + ax)*sx 

    #time dependent
    Hcot = v*((w/g) + ((mu**2)/(g*w)))*(ax.dag() + ax)*sx 
    
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
    nS_transition_txt = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\data_new\\nS_transition_d={d}_mu={mu}_g={g}'
    np.savetxt(nS_transition_txt, np.c_[eigvl, ns_trans])
