import numpy as np
import matplotlib.pyplot as plt
from mm_parameters import*
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm

# This code read data from files with energies and transition rates from the |nS> state for each quasienergy operator size, given by the number of blocks.

# Read data from each file and store in lists or arrays
#dx = [1]
Nb = 49
data = []

# Load data from files
for k in range(0, 351):
    file_path = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\data\\nS_transition_dx={k}_N=50_Nb=49_gx=4'
    with open(file_path, 'r') as file:
        lines = file.readlines()
        E_values = []
        nS_values = []
        for line in lines:
            x, y = map(float, line.split())
            E_values.append(x - 250)
            nS_values.append(y)
        data.append((E_values, [k] * len(E_values), nS_values))

# Create structured data for 3D plotting
E_values_all = []
dx_values_all = []
nS_values_all = []

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24

for E_values, dx_values, nS_values in data:
    E_values_all.extend(E_values)
    dx_values_all.extend(dx_values)
    nS_values_all.extend(nS_values)


# Sort data based on the z-coordinate (nS_values_all) in descending order
sorted_indices = np.argsort(nS_values_all)[::1]
E_values_sorted = np.array(E_values_all)[sorted_indices]
dx_values_sorted = np.array(dx_values_all)[sorted_indices]
nS_values_sorted = np.array(nS_values_all)[sorted_indices]

# Define your tolerance
tolerance = 2e-2  # Adjust this value according to your preference

# Clip values below the tolerance
nS_values_clipped = np.clip(nS_values_sorted, tolerance, None)

# Plot with clipped z-values
fig2 = plt.figure(figsize=(8, 5))
ax2 = fig2.add_subplot(111)
nS_trans = ax2.scatter(dx_values_sorted/hol, E_values_sorted, c=nS_values_clipped, norm=LogNorm(), cmap='binary', s=1, rasterized=True)

# Set limits for the y-axis (energy values) in the projection plot
ax2.set_xlim(0, 35)
ax2.set_ylim(-50200, -49800)

# Set labels and title for the projection plot
ax2.set_xlabel(r'$\kappa_{x}$', fontsize=24)
ax2.set_ylabel(r'$\epsilon_{\alpha m}$ (units of $\hbar \omega$)', fontsize=24)

nS_trans_bar = plt.colorbar(nS_trans)
nS_trans_bar.set_label(r'$\left|\langle nS|u_{nm} \rangle\right|^2$', fontsize=24)
plt.tight_layout()

#fig_save = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\figures\\quasi_spectrum_gx4.pdf'
#plt.savefig(fig_save, format='pdf')
plt.show()