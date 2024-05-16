import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
import time
from mm_operators import*
#from mm_operators1d import*
from mm_parameters import*

start_time = time.time()

# calculating the Floquet quasi-energies and Floquet modes associated with the system
dx_vec = np.linspace(0, 1.6, 40)

hs_dim = Nx*Ny*2

q_energies = np.zeros((len(dx_vec), hs_dim))
ns_trans = np.zeros((len(dx_vec), hs_dim))
T = 2 * np.pi / v 

for idx, dx in enumerate(dx_vec):

    #shorthand notations
    kx, ky = dx/hol, dy/hol
    
    ############### internal motion #################
    Hin = E*sz + 0.25*rho*(kx*sx + ky*sy)

    ############### external motion #################
    Hex0 = gx*ax.dag()*ax + gy*ay.dag()*ay + (kx*ux*np.sqrt(gx))*(ax.dag() + ax) + (ky*uy*np.sqrt(gy))*(ay.dag() + ay)

    #time dependent on sine and cosine
    Hex1t = -1j*((zx*gx)*(ax.dag()*ax.dag() - ax*ax) - zy*gy*(ay.dag()*ay.dag() - ay*ay))
    Hex2t = 0.5*(zx**2)*gx*(ax.dag()*ax.dag() + 2*ax.dag()*ax + ax*ax) + 0.5*(zy**2)*gy*(ay.dag()*ay.dag() + 2*ay.dag()*ay + ay*ay)

    ############### coupling dynamics ###############
    Hco0 = -rho*(ux*np.sqrt(gx)*(ax.dag() + ax)*sx + uy*np.sqrt(gy)*(ay.dag() + ay)*sy) 

    #time dependent
    Hcot = rho*xi*((zx*np.sqrt(gx))*(ax.dag() + ax)*sx - (zy*np.sqrt(gy))*(ay.dag() + ay)*sy)

    #functions sine and cosine utilized in the code
    sinvt, cosvt, cos2vt = 'np.sin({}*t)', 'np.cos({}*t)', 'np.cos({}*2*t)' 

    sinvt, cosvt, cos2vt = sinvt.format(v), cosvt.format(v), cos2vt.format(v) 

    #generating time dependent Hamiltonians
    H0 = Hin + Hex0 + Hco0

    H = [H0, [Hex1t, sinvt], [Hex2t, cos2vt], [Hcot, cosvt]]

    #state |nS> transition 
    nS_state = qt.tensor(qt.basis(2, 0), qt.destroy(Nx), qt.qeye(Ny)) 

    opts = qt.Options(nsteps=500000)
    f_modes, f_energies = qt.floquet_modes(H, T, options=opts)
    q_energies[idx, :] = f_energies
    ns_trans[idx, :] = nS_state.dag()*f_modes

################# Plotting of the transitions to the state nS #######################
plt.figure(figsize=(5,5))

plt.rcParams['mathtext.fontset'] = 'stix'
# Set the font family and size
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24

#plt.xlabel(r'$\delta_x$', fontsize=24)
plt.ylabel(r'$\left|<nS|\epsilon_{nm}>\right|$', fontsize=24)

plt.plot(dx_vec, ns_trans[:, 10], color='black', marker='.', linestyle='',
            markersize=1)

plt.show()

# Plotting
plt.figure(figsize=(5.5, 12))  # Adjust the figure size as needed

plt.rcParams['mathtext.fontset'] = 'stix'
# Set the font family and size
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24

plt.xlabel(r'$\delta_x$', fontsize=24)
plt.ylabel('Quasienergies', fontsize=24)

        
# Define the frequency zones and their corresponding colors
frequency_zones = {
    (-1.5, -0.5): '#5ab4ac',  # Define the frequency zones and their colors here
    (-0.5, 0.5): '#c7eae5',
    (0.5, 1.5): '#5ab4ac'
}

# Setting plot limits
plt.xlim(0, 1.6) # setting y-axis limits from 0 to 12

# Plotting
for m in range(-1, 2):
    for k in range(0, hs_dim):
        plt.plot(dx_vec, q_energies[:, k] + m*(v/w), color='black', marker='.', linestyle='',
            markersize=1)

# Add background patches for each frequency zone
for zone, color in frequency_zones.items():
    plt.axhspan(zone[0] * v / w, zone[1] * v / w, facecolor=color, alpha=0.3)

# Set y-ticks and labels
y_ticks = np.arange(-1, 2) * v / w
y_labels = ['$\\epsilon-\\hbar\\nu$', '$\epsilon$', '$\\epsilon+\\hbar\\nu$']  # You can modify these labels as needed

plt.yticks(y_ticks, y_labels)

plt.tight_layout(pad=1.0)  # Adjust padding as needed
#plt.savefig('floquet_quasienergies_test.pdf', format='pdf')
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))
 