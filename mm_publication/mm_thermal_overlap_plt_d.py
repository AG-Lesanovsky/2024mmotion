import numpy as np
import matplotlib.pyplot as plt

# Define parameters
eta = 10
n_max = 25  # Set the maximum n value
x_list = range(0, 201)
omega = 1.0  # Frequency (set your actual value)
T = 5  # Temperature (set your actual value)

# Set up the figure
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24

fig = plt.figure(figsize=(7, 5))
ax = fig.add_subplot(111)

# Lists to store all data
E_values_all, dx_values_all, nS_values_all, boltzmann_factors = [], [], [], []

# Compute the partition function Z
Z = sum(np.exp(-n * omega / T) for n in range(n_max + 1))

# Loop over n values
for n in range(n_max + 1):
    boltzmann_factor = np.exp(-n * omega / T) / Z  # Boltzmann weight for each n
    for x in x_list:
        d = x / 2
        file_path = f'C:/Users/Wilson/OneDrive/Dokumente/2024mmotion/1d_model/finite_temp/data_finite/overlap_d={d}_eta={eta}_n={n}'
        try:
            with open(file_path, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    E, nS = map(float, line.split())
                    E_values_all.append(E - 2*n)
                    dx_values_all.append(d)
                    nS_values_all.append(nS * boltzmann_factor)  # Apply Boltzmann weight
        except FileNotFoundError:
            print(f"File not found for n = {n}, d = {d}")

# Convert to numpy arrays
E_values_all = np.array(E_values_all)
dx_values_all = np.array(dx_values_all)
nS_values_all = np.array(nS_values_all)

# Sort data by overlap magnitude
sorted_indices = np.argsort(nS_values_all)
E_values_sorted = E_values_all[sorted_indices]
dx_values_sorted = dx_values_all[sorted_indices]
nS_values_sorted = nS_values_all[sorted_indices]

# Define tolerance
tolerance = 5e-5

# Clip values below the tolerance
nS_values_clipped = np.clip(nS_values_sorted, tolerance, None)

# Filter points where nS_values_clipped > 0.02
filtered_indices = np.where(nS_values_clipped > 0.02)[0]

# Extract the filtered data
E_values_filtered = E_values_sorted[filtered_indices]
dx_values_filtered = dx_values_sorted[filtered_indices]
nS_values_filtered = nS_values_clipped[filtered_indices]

# Scatter plot
nS_trans = ax.scatter(dx_values_filtered, E_values_filtered, c=nS_values_filtered, cmap='bone_r', s=5, rasterized=True, vmin=0, vmax=0.5, label = r'$\langle N \rangle = 0.1, T = 0.42$')
plt.yticks([-30010, -30000, -29990], [r'$-10$', r'$0$', r'$10$'])
plt.xticks([0, 50, 100], [r'$0$', r'$50$', r'$100$'])

ax.set_xlim(0, 100)
ax.set_ylim(-30010, -29990)
ax.set_xlabel(r'$\delta$', fontsize=24)
ax.set_ylabel(r'$E$', fontsize=24)

plt.text(10, -30000, r'(a) $\zeta = 3.5$', fontsize = 24)

ax.legend(loc = 4, frameon=False)

# Colorbar
nS_trans_bar = plt.colorbar(nS_trans)
nS_trans_bar.set_label(r'$P_{E}$', fontsize=24)

plt.tight_layout()
# fig_save = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\finite_temp\\figures_finite\\thermal_z={eta}_T={T}.pdf'
# plt.savefig(fig_save, format='pdf')
plt.show()






