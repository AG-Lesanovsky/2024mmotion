import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
import time
from mm_operators import*
#from mm_operators1d import*
from mm_parameters import*

start_time = time.time()

# calculating the Floquet quasi-energies and Floquet modes associated with the system
hs_dim = Nx*Ny*2

ns_trans = np.zeros(hs_dim)
T = 2 * np.pi / v 

xspan = np.linspace(0, 100, hs_dim)

#state |nS> transition 
nS_state = qt.tensor(qt.basis(2, 0), qt.basis(Nx, 0), qt.basis(Nx, 0)) 
#projector unto |nS>
PnS = nS_state*nS_state.dag()

opts = qt.Options(nsteps=500000)
f_modes, f_energies = qt.floquet_modes(H, T, options=opts)

for k in range(hs_dim):
    ns_trans[k] = qt.expect(PnS, f_modes[k])

#print(nS_state)
#print(f_energies)
#print(nS_state.dag()*f_modes[0])

print("--- %s seconds ---" % (time.time() - start_time))

################## Plotting of the transitions to the state nS #######################
plt.figure(figsize=(7,5))

plt.rcParams['mathtext.fontset'] = 'stix'
# Set the font family and size
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24

plt.ylim(0, 1)

plt.xlabel(r'$\epsilon_{nm} (\omega)$', fontsize=24)
plt.ylabel(r'$\left| \langle nS|u_{nm} \rangle \right|^2$', fontsize=24)

#plt.plot(f_energies, ns_trans, color='black', marker='.', linestyle='',
#            markersize=4)

plt.vlines(f_energies, 0, ns_trans, color='#018571', linestyle='solid', linewidth = 1.5)

plt.tight_layout()

#saving plots
nS_transition_pdf = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\figures\\nS_transition_dx={dx}.pdf'
plt.savefig(nS_transition_pdf, format='pdf')

#saving data
nS_transition_txt = f'C:\\Users\\Wilson\\OneDrive\\Dokumente\\2024mmotion\\data\\nS_transition_dx={dx}'
np.savetxt(nS_transition_txt, np.c_[f_energies, ns_trans])

plt.show()