import numpy as np
import matplotlib.pyplot as plt
#from mm_parameters import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize, LogNorm

#The present code contains the tools used for plotting figures of the colourmaps presented in Fig. 2(a-c), for variable \delta 
#The entries are the quasienergies and respective overlap magnitudes for different frequency \delta

# Read data from each file and store in lists or arrays
data = []
x_values = np.linspace(0, 300, 1000)

eta = 3.5
n_blocks = 21
n = 1

x_list = range(0, 201)

#for x in range(0, 201):
for x in x_list:
    # d = x/2
    d = x/2
    file_path = f'C:/Users/Wilson/OneDrive/Dokumente/2024mmotion/1d_model/finite_temp/data_finite/overlap_d={d}_eta={eta}_n={n}'
    with open(file_path, 'r') as file:
        lines = file.readlines()
        E_values = []
        nS_values = []
        for line in lines:
            x, y = map(float, line.split())
            #E_values.append(x-71*eta)
            E_values.append(x - 2*n)
            nS_values.append(y)
            
        data.append((E_values, [d] * len(E_values), nS_values))

 

# Calculate the energy function values
d_list = np.linspace(0, 100, 1000)

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

# Now filter points where nS_values_clipped > 0.5
filtered_indices = np.where(nS_values_clipped > 0.02)[0]

# Extract the filtered data
E_values_filtered = E_values_sorted[filtered_indices]
dx_values_filtered = dx_values_sorted[filtered_indices]
nS_values_filtered = nS_values_clipped[filtered_indices]


fig = plt.figure(figsize=(7, 5))
ax = fig.add_subplot(111)
nS_trans = ax.scatter(dx_values_sorted, E_values_sorted, c=nS_values_clipped, cmap='bone_r', s=5, rasterized=True,vmin=0, vmax=1, label = r'$\eta = 10$')

ax.set_xlim(0, 100)
ax.set_ylim(-30010, -29990)
# ax.set_ylim(29960, 30000)


ax.set_xlabel(r'$\delta$', fontsize=24)
ax.set_ylabel(r'$E$')

ax.legend(loc = 4, frameon=False)

nS_trans_bar = plt.colorbar(nS_trans)
#nS_trans_bar.set_label(r'overlap $\left|\langle\langle \downarrow|u_{\alpha, m}\rangle \rangle\right|^2$', fontsize=24)
# plt.yticks([-30000.5, -30000, -29999.5], [r'$-0.5$', r'$0$', r'$0.5$'])
nS_trans_bar.set_label(r'$P_{E}$', fontsize=24)
#plt.yticks([-10, 0, 10], [r'$-10$', r'$0$', r'$10$'])
#plt.xticks([0, 50], [r'$0$', r'$50$'])


# plt.text(10, -30000, '(a)', fontsize = 24)
#plt.text(90, 8, '(c)', fontsize = 24)
plt.tight_layout()

# fig_save = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\finite_temp\\figures_finite\\N={n}_eta={eta}.pdf'
# plt.savefig(fig_save, format='pdf')
plt.show()
