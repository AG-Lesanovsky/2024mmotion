import numpy as np
import matplotlib.pyplot as plt
#from mm_parameters import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize, LogNorm

#The present code contains the tools used for plotting figures of the colourmaps presented in Fig. 2(d-f), for variable \eta 
#The sntries are the quasienergies and respective overlap magnitudes for different frequency \eta

#Read data from each file and store in lists or arrays
data = []
x_values = np.linspace(0, 300, 1000)

d = 0
n_blocks = 26
n = 0



for x in range(0, 201):
    eta = 3.5 + x/25
    # eta = 8 + 0.2*x
    #file_path = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\data\\smnS_transition_d={d}_eta={eta}'
    file_path = f'C:/Users/wsant/OneDrive/Dokumente/2024mmotion/1d_model/finite_temp/data_finite/overlap_d={d}_eta={eta}_n={n}'
    with open(file_path, 'r') as file:
        lines = file.readlines()
        E_values = []
        nS_values = []
        for line in lines:
            x, y = map(float, line.split())
            E_values.append(x)
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


fig = plt.figure(figsize=(7, 5))
#fig = plt.figure(figsize=(4, 3))
ax = fig.add_subplot(111)
nS_trans = ax.scatter(dx_values_sorted, E_values_sorted, c=nS_values_clipped, cmap='bone_r', s=5, rasterized=True, vmin=0, vmax=1, label = r'$\delta = 100$')
# ax.plot(eta_list, energy_2_vareta, linestyle='dashed', color='#b2182b', linewidth=3)

ax.set_xlim(3.5, 11.5)
ax.set_ylim(-30010, -29990)
#ax.set_ylim(-10, 10)

ax.set_xlabel(r'$\eta$', fontsize=24)
ax.set_ylabel(r'$E$')

ax.legend(loc = 4, frameon=False)

nS_trans_bar = plt.colorbar(nS_trans)
nS_trans_bar.set_label(r'$P_{E}$', fontsize=24)
#plt.yticks([-29990, -30000, -30010], [r'$10$', r'$0$', r'$-10$'])
#plt.yticks([10, 0, -10], [r'$10$', r'$0$', r'$-10$'])

#plt.text(9.8, -29999, '(d)', fontsize = 24)
#plt.text(9, 8, '(f)', fontsize = 24)
plt.tight_layout()

#fig_save = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\figures\\quasi_spectrum_d={d}.pdf'
# fig_save = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\figures\\ground_state_delta={d}.pdf'
# plt.savefig(fig_save, format='pdf')
plt.show()
