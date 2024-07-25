import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize, LogNorm
from mm_overlap_pert import *

# This code reads data from files with energies and transition rates from the |nS> state for each quasienergy operator size, given by the number of blocks.

# Read data from each file and store in lists or arrays
Nb = 49
data = []
x_values = np.linspace(0, 300, 1000)

eta = 3.5
#TRY TO REMEMBER THAT FOR eta = 3.5 YOU HAVE USED N = 36

for x in range(0, 201):
    d = x/2
    file_path = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\data_new\\smnS_transition_d={d}_eta={eta}'
    with open(file_path, 'r') as file:
        lines = file.readlines()
        E_values = []
        nS_values = []
        for line in lines:
            x, y = map(float, line.split())
            E_values.append(x - 21*eta)
            nS_values.append(y)
        data.append((E_values, [d] * len(E_values), nS_values))

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

# Define tolerance
tolerance = 5e-5

ovlp_m00_d = 1/norm(delta_list, eta)

# Clip values below the tolerance
nS_values_clipped = np.clip(nS_values_sorted, tolerance, None)

fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(111)
ax.set_xlabel(r'misalignment length $\delta$', fontsize=24)
ax.set_ylabel(r'overlap $\left|\langle \langle \downarrow|\downarrow^{(1)} \rangle \rangle\right|^{2}$')
ax.plot(dx_values_sorted, nS_values_clipped, color='#2166ac', linewidth=3)
ax.plot(delta_list, norm_m00_dx**2, linestyle='dashed', color='#b2182b', linewidth=3, label = r'$\eta = 3.5$')

ax.legend(loc = 4, frameon=False)

ax.set_xlim(0, 100)
ax.set_ylim(0, 1)


plt.text(4, 0.1, '(c)', fontsize = 24)
plt.tight_layout()

fig_save = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\figures\\pert_overlap_eta={eta}.pdf'
plt.savefig(fig_save, format='pdf')
plt.show()
