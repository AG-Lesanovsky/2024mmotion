import numpy as np
import matplotlib.pyplot as plt
from mm_parameters import*
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
from matplotlib.colors import LogNorm

# This code read data from files with energies and transition rates from the |nS> state for each quasienergy operator size, given by the number of blocks.

# Read data from each file and store in lists or arrays
data = []

# Load data from files
for k in range(0, 351):
    file_path = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\data\\nS_transition_dx={k}_N=50_Nb=49'
    with open(file_path, 'r') as file:
        lines = file.readlines()
        E_values = []
        nS_values = []
        for line in lines:
            x, y = map(float, line.split())
            E_values.append(x)
            nS_values.append(y)
        data.append((E_values - np.rint(0.5*(49 - 1))*xi, [k] * len(E_values), nS_values))

# Create structured data for 3D plotting
E_values_all = []
dx_values_all = []
nS_values_all = []

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24

# Filter data based on x-axis limits
#filtered_data = [(E_values, dx_values, nS_values) for E_values, dx_values, nS_values in data if min(E_values) >= -50010 and max(E_values) <= -49950]

for E_values, dx_values, nS_values in data:
    E_values_all.extend(E_values)
    dx_values_all.extend(dx_values)
    nS_values_all.extend(nS_values)

# Plot 3D scatter plot
fig1 = plt.figure(figsize=(10, 8))
# Set background color
ax1 = fig1.add_subplot(111, projection='3d')
ax1.scatter(E_values_all, dx_values_all/hol, nS_values_all, c=nS_values_all, cmap='Greys')
# First remove fill
ax1.xaxis.pane.fill = True
ax1.yaxis.pane.fill = True
ax1.zaxis.pane.fill = True

# Now set color to white (or whatever is "invisible")
#ax1.xaxis.pane.set_edgecolor('black')
#ax1.yaxis.pane.set_edgecolor('black')
#ax1.zaxis.pane.set_edgecolor('black')

# Bonus: To get rid of the grid as well:
#ax1.grid(False)

# Set limits for the x-axis (energy values)
ax1.set_xlim(-50010, -49950)

# Set labels and title
ax1.set_xlabel(r'$\epsilon_{\alpha m}(\hbar \omega)$', fontsize=24)
ax1.set_ylabel(r'$\kappa_{x}$ (nm)', fontsize=24)
ax1.set_zlabel(r'$\left|\langle nS|u_{\alpha m} \rangle\right|^2$', fontsize=24)

plt.tight_layout()

plt.show()

