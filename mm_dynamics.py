#code for running the micromotion using Floquet dynamics
#all the variables are obtained and defined as in Joe's text
#troughout the text we use LaTeX notation for Greek letters & symbols  
#everything is in units of hbar times trap frequency \omega
import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
import time
from mm_operators import*
from mm_parameters import*

start_time = time.time()

#the domain of integration
tmax = 50
npts = 1000
tspan = np.linspace(0, tmax, npts)

#initial state
phi0 = qt.tensor(qt.basis(4, 0), qt.basis(Nx, 0), qt.basis(Ny, 0), qt.basis(Nz, 0))

result = qt.sesolve(H, phi0, tspan, [ax.dag()*ax, ay.dag()*ay, az.dag()*az.dag()])

#expected values for each operators
ax_exp = result.expect[0]
ay_exp = result.expect[1]
az_exp = result.expect[2]

print("--- %s seconds ---" % (time.time() - start_time))

plt.plot(tspan, ax_exp)
plt.plot(tspan, ay_exp)
plt.plot(tspan, az_exp)
plt.show()