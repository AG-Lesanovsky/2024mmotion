import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize, LogNorm
from mm_overlap_pert import *

#The present code contains the tools used for plotting figures of the curves presented in Fig. 3(d-f), for variable \eta 
#The sntries are the quasienergies and respective overlap magnitudes for different frequency \eta
#We also use the functions from the script mm_overlap_pert.py, containing the matrix elements and respective perturbative energies

# Read data from each file and store in lists or arrays
Nb = 49
data = []
x_values = np.linspace(0, 300, 1000)

d = 100

for x in range(0, 201):
    eta = 3.5 + x/25
    file_path = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\data\\smnS_transition_d={d}_eta={eta}'
    with open(file_path, 'r') as file:
        lines = file.readlines()
        E_values = []
        nS_values = []
        for line in lines:
            x, y = map(float, line.split())
            E_values.append(x - 36*eta)
            nS_values.append(y)
        data.append((E_values, [eta] * len(E_values), nS_values))

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

# Clip values below the tolerance
nS_values_clipped = np.clip(nS_values_sorted, tolerance, None)

fig = plt.figure(figsize=(6, 5))
ax = fig.add_subplot(111)
ax.set_xlabel(r'dimensionless frequency $\eta$', fontsize=24)
ax.set_ylabel(r'overlap $\langle \langle \downarrow|u_{\downarrow} \rangle \rangle$')
ax.plot(dx_values_sorted, nS_values_clipped,  color='#2166ac', linewidth=3)
ax.plot(eta_list, norm_m00_eta**2, linestyle='dashed',  color='#b2182b', linewidth=3, label = r'$\langle \langle \downarrow|\downarrow^{(2)} \rangle \rangle$')

ax.legend(loc = 4, frameon=False)

ax.set_xlim(3.5, 11.5)
ax.set_ylim(0, 1)

plt.text(4, 0.1, '(f)', fontsize = 24)
plt.tight_layout()

fig_save = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\figures\\pert_overlap_delta={delta}.pdf'
plt.savefig(fig_save, format='pdf')
plt.show()

