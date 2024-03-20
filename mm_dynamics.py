#script to run the micromotion using Floquet dynamics
#all the variables are obtained and defined as in Joe's text
#troughout the text we use LaTeX notation for Greek letters & symbols  
#everything is in units of hbar times trap frequency \omega
import qutip as qt
import numpy as np
import matplotlib.pyplot as plt
import time
#from mm_operators import*
from mm_operators1d import*
from mm_parameters import*

start_time = time.time()

#the domain of integration
tmax = 40
npts = 1000000
tspan = np.linspace(0, tmax, npts)

#initial state
#phi0 = qt.tensor(qt.basis(4, 0), qt.basis(Nx, 0), qt.basis(Ny, 0), qt.basis(Nz, 0))
phi0 = qt.tensor(qt.basis(2, 0), qt.basis(Nx, 0))

#result = qt.sesolve(H, phi0, tspan, [ax.dag()*ax, ay.dag()*ay, az.dag()*az.dag()])
result = qt.sesolve(H, phi0, tspan, [ax.dag()*ax, S1*S1d, P1*P1d])

#expected values for each operators
ax2_exp = result.expect[0]
s1_exp = result.expect[1]
p1_exp = result.expect[2]
#ay_exp = result.expect[1]
#az_exp = result.expect[2]

print("--- %s seconds ---" % (time.time() - start_time))

#plt.plot(tspan, ax2_exp)
plt.plot(tspan, (1- s1_exp)*1e5)
#plt.plot(tspan, ay_exp)
#plt.plot(tspan, az_exp)
plt.show()

plt.plot(tspan, p1_exp*1e5)
#plt.plot(tspan, ay_exp)
#plt.plot(tspan, az_exp)
plt.show()

plt.plot(tspan, ax2_exp)
#plt.plot(tspan, ay_exp)
#plt.plot(tspan, az_exp)
plt.show()