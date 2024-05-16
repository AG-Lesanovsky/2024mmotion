import numpy as np
import matplotlib.pyplot as plt
from mm_parameters import*
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize

# This code read data from files with energies and transition rates from the |nS> state for each quasienergy operator size, given by the number of blocks.

# Read data from each file and store in lists or arrays
dx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
Nb = 51
data = []

# Load data from files
for k in dx:
    file_path = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\2d_model\\data\\nS_transition_dx={k}_N=50_Nb=51'
    with open(file_path, 'r') as file:
        lines = file.readlines()
        E_values = []
        nS_values = []
        for line in lines:
            x, y = map(float, line.split())
            E_values.append(x)
            nS_values.append(y)
        data.append((E_values - np.rint(0.5*(Nb - 1))*xi, [k] * len(E_values), nS_values))

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

# Plot 3D scatter plot
fig1 = plt.figure(figsize=(8, 8))
# Set background color
ax1 = fig1.add_subplot(111, projection='3d')
ax1.scatter(E_values_all, dx_values_all, nS_values_all, c=nS_values_all, cmap='Greys')
# First remove fill
ax1.xaxis.pane.fill = True
ax1.yaxis.pane.fill = True
ax1.zaxis.pane.fill = True

# Now set color to white (or whatever is "invisible")
ax1.xaxis.pane.set_edgecolor('black')
ax1.yaxis.pane.set_edgecolor('black')
ax1.zaxis.pane.set_edgecolor('black')

# Bonus: To get rid of the grid as well:
ax1.grid(False)

# Set limits for the x-axis (energy values)
ax1.set_xlim(-50000, -49980)

# Set labels and title
ax1.set_xlabel(r'$\epsilon_{\alpha m}(\omega)$', fontsize=24)
ax1.set_ylabel(r'$\delta_{x}$', fontsize=24)
ax1.set_zlabel(r'$|\langle nS|u_{\alpha m} \rangle|^2$', fontsize=24)

# Sort data based on the z-coordinate (nS_values_all) in descending order
sorted_indices = np.argsort(nS_values_all)[::1]
E_values_sorted = np.array(E_values_all)[sorted_indices]
dx_values_sorted = np.array(dx_values_all)[sorted_indices]
nS_values_sorted = np.array(nS_values_all)[sorted_indices]

# Define the maximum value for the z-axis
z_max = 0.6  # Adjust this value according to your preference

# Replace z-values equal to or greater than z_max with z_max
nS_values_clipped = np.clip(nS_values_sorted, None, z_max)

# Plot projection onto x-y plane (cut from above) with clipped z-values
fig2 = plt.figure(figsize=(5, 5))
ax2 = fig2.add_subplot(111)
ax2.scatter(dx_values_sorted, E_values_sorted, c=nS_values_clipped, cmap='Greys')

# Set logarithmic scale along the x-axis
ax2.set_xscale('log')

# Set limits for the y-axis (energy values) in the projection plot
ax2.set_ylim(-50000, -49980)

# Set labels and title for the projection plot
ax2.set_xlabel(r'$\delta_{x}$', fontsize=24)
ax2.set_ylabel(r'$\epsilon_{\alpha m}(\omega)$', fontsize=24)

plt.tight_layout()

plt.show()
