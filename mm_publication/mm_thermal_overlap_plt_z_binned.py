import numpy as np
import matplotlib.pyplot as plt

# Define parameters
d = 100
n_max = 25  # Set the maximum n value
x_list = range(0, 201)
omega = 1.0  # Frequency (set your actual value)
T = 0.42 # Temperature (set your actual value)

# Set up the figure
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24

fig = plt.figure(figsize=(7, 5))
ax = fig.add_subplot(111)

# Lists to store all data
data_dict = {}

# Compute the partition function Z
Z = sum(np.exp(-n * omega / T) for n in range(n_max + 1))

# Loop over n values
for n in range(n_max + 1):
    boltzmann_factor = np.exp(-n * omega / T) / Z  # Boltzmann weight for each n
    for x in x_list:
        eta = 3.5 + x/25
        file_path = f'C:/Users/Wilson/OneDrive/Dokumente/2024mmotion/1d_model/finite_temp/data_finite/overlap_d={d}_eta={eta}_n={n}'
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    E, nS = map(float, line.split())
                    E_adj = E - 2 * n
                    if -30010 <= E_adj <= -29990:  # Only consider energies in the desired range
                        if eta not in data_dict:
                            data_dict[eta] = []
                        data_dict[eta].append((E_adj, nS * boltzmann_factor))  # Store energy and weighted overlap
        except FileNotFoundError:
            print(f"File not found for n = {n}, eta = {eta}")

# Define energy bins within the desired range
bin_width = 1.01  # Adjust as needed
E_min, E_max = -30011, -29989  # Restrict binning to this range
E_bins = np.arange(E_min, E_max + bin_width, bin_width)

# Prepare data for plotting
binned_dx = []
binned_E = []
binned_nS = []

for eta, values in data_dict.items():
    E_values, nS_values = zip(*values)  # Unpack energies and overlaps
    E_values = np.array(E_values)
    nS_values = np.array(nS_values)
    
    binned_nS_values = np.zeros(len(E_bins) - 1)
    binned_center_of_mass = np.zeros(len(E_bins) - 1)
    
    for i in range(len(E_bins) - 1):
        in_bin = (E_values >= E_bins[i]) & (E_values < E_bins[i + 1])
        # Weighted sum for center of mass calculation
        if np.sum(nS_values[in_bin]) > 0:  # Avoid division by zero
            binned_center_of_mass[i] = np.sum(E_values[in_bin] * nS_values[in_bin]) / np.sum(nS_values[in_bin])
        binned_nS_values[i] = np.sum(nS_values[in_bin])
    
    # Store nonzero binned data
    for i, binned_val in enumerate(binned_nS_values):
        if binned_val > 0:
            binned_dx.append(eta)
            binned_E.append(binned_center_of_mass[i])  # Use center of mass for the bin
            binned_nS.append(binned_val)

# Convert lists to arrays
binned_dx = np.array(binned_dx)
binned_E = np.array(binned_E)
binned_nS = np.array(binned_nS)

# Plot the binned data
nS_trans = ax.scatter(binned_dx, binned_E, c=binned_nS, cmap='bone_r', s=5, rasterized=True, vmin=0, vmax=1, label=r'$\langle N \rangle = 0.1, T = 0.42$')
plt.yticks([-30010, -30000, -29990], [r'$-10$', r'$0$', r'$10$'])
plt.xticks([5, 10], [r'$5$', r'$10$'])

ax.set_xlim(3.5, 11.5)
ax.set_ylim(-30010, -29990)
ax.set_xlabel(r'$\zeta$', fontsize=24)
ax.set_ylabel(r'$E + \varepsilon$', fontsize=24)

plt.text(5, -29995, r'(b) $\delta = 100$', fontsize=24)
ax.legend(loc=4, frameon=False)

# Colorbar
nS_trans_bar = plt.colorbar(nS_trans)
nS_trans_bar.set_label(r'$P_{E}(T)$', fontsize=24)

plt.tight_layout()
fig_save = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\finite_temp\\figures_finite\\thermal_d={d}_T={T}.pdf'
plt.savefig(fig_save, format='pdf')
plt.show()
