import matplotlib.pyplot as plt
import numpy as np

data_paths = [
    'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\data\\op_norm_dxvar'
]

# Initialize empty lists to store eigenvalues and nS_trans values
dx_list = []
op_norm_values = []

# Loop through each data path
for path in data_paths:
    # Read data from the file
    with open(path, 'r') as file:
        for line in file:
            # Split the line into eigenvalue and nS_trans value (assuming they are space-separated)
            values = line.split()
            dx_list.append(float(values[0])/10)  # Assuming eigenvalue is the first value
            op_norm_values.append(float(values[1]))  # Assuming nS_trans is the second value

    
plt.figure(figsize=(6,5))

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24

plt.xlim(0, 4e4)
plt.ylim(0, 1)
plt.xlabel(r'$\kappa_{x}$', fontsize=24)
plt.ylabel(r'$\left\|\mathbb{V}\right\|/\left\|\mathbb{H}_\mathrm{F}\right\|$', fontsize=24)

# Plot the vertical lines
plt.plot(dx_list, op_norm_values, color='#b2182b', linestyle='solid', linewidth=1.5, marker='^', markersize=10)
plt.tight_layout()
op_norm_dxvar = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\figures\\op_norm_dxvar.pdf'
plt.savefig(op_norm_dxvar, format='pdf')

#plt.xscale('log')


plt.show()

#saving plots
