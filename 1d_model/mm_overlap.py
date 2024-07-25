import matplotlib.pyplot as plt
import numpy as np

d = 50.0
eta = 10

data_paths = [
    f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\data\\smnS_transition_d={d}_eta={eta}'
]

# Initialize empty lists to store eigenvalues and nS_trans values
qe_list = []
nS_values = []

# Loop through each data path
for path in data_paths:
    # Read data from the file
    with open(path, 'r') as file:
        for line in file:
            # Split the line into eigenvalue and nS_trans value (assuming they are space-separated)
            values = line.split()
            qe_list.append(float(values[0]) - 21*eta)  # Assuming eigenvalue is the first value
            nS_values.append(float(values[1]))  # Assuming nS_trans is the second value

    
plt.figure(figsize=(3,2.5))

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24

# Plot the vertical lines
plt.vlines(qe_list, 0, nS_values, color='#b2182b', linestyle='solid', linewidth = 5)
plt.tight_layout()

# Set axis limits
plt.xlim([-30010, -29990])  # Set the limits for the x-axis
plt.ylim([0.01, 1])
plt.yscale('log')
plt.xticks(labels=None)

nS_tran = f'C:\\Users\\wsant\\OneDrive\\Dokumente\\2024mmotion\\1d_model\\figures\\nS_trans_d={d}_eta={eta}.pdf'
plt.savefig(nS_tran, format='pdf')

plt.show()

#saving plots
