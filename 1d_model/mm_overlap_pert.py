import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize, LogNorm

#In this code we calculate curves with results obtained via perturbation theory. For clarity, we define each part of the code and structuralize it accordingly
E = 30000
delta = 50
eta = 10
eta_list = np.linspace(0, 15, 1000)
delta_list = np.linspace(0, 100, 1000)

#energies at zeroth order with respect to index, frequency 
def energy_0(s, n, m, eta):
    return s*E + n + m*eta

#energy of the state \downarrow, 0, 0 at first order with respect to frequency and displacement
def energy_2(delta, eta):
    a0 = (((9*delta)/32)**2)*(1/(2*E)) + (1/(16*(2*E + 1))) 
    a1 = (((3*delta)/64)**2)*(2/(1-eta**2)) + (((3*np.sqrt(2))/8)**2)*(4/(4-eta**2)) + (((3*eta*delta)/16)**2)*((4*E)/(4*E**2 - eta**2)) + ((9*eta**2)/4)*((4*E + 2)/((2*E + 1)**2 - eta*2))
    a2 = (((9*delta)/256)**2)*(2/(1 - 4*eta**2)) + (((9*np.sqrt(2))/64)**2)*(1/(1 - eta**2))
    return - E - a0 - a1 - a2

#true normalized state
def norm_m00(delta, eta):
    #matrix elements 
    mep00 = -(9*delta)/32
    mep10 = 1/4
    mem01 = -3/8
    mem11 = (3*delta)/64
    mem21 = -(3*np.sqrt(2))/8
    mep01 = (3*eta*delta)/16
    mep11 = -(3*eta)/2 
    mem02 = - 9/64
    mem12 = -(9*delta)/256
    mem22 = -(9*np.sqrt(2))/64
    
    #energies associated to each transition
    Ep00 = 2*E
    Ep10 = 2*E + 1

    Em01 = eta
    Em11 = 1 + eta
    Em21 = 2 + eta
    Em01m = - eta
    Em11m = 1 - eta
    Em21m = 2 - eta
    Ep01 = 2*E + eta
    Ep11 = 2*E + 1 + eta
    Ep01m = 2*E - eta
    Ep11m = 2*E + 1 - eta 
    Em02 = 2*eta
    Em12 = 1 + 2*eta
    Em22 = 2 + 2*eta
    Em02m = -2*eta
    Em12m = 1 - 2*eta
    Em22m = 2 - 2*eta

    #sum
    norm_c = (mep00/Ep00)**2 + (mep10/Ep10)**2 + (mem01/Em01)**2 + (mem01/Em01m)**2
    norm_c = norm_c + (mem11/Em11)**2 + (mem11/Em11m)**2 + (mem21/Em21)**2 + (mem21/Em21m)**2
    norm_c = norm_c + (mep01/Ep01)**2 + (mep01/Ep01m)**2 + (mep11/Ep11)**2 + (mep11/Ep11m)**2 + (mem02/Em02)**2 + (mem02/Em02m)**2
    norm_c = norm_c + (mem12/Em12)**2 + (mem12/Em12m)**2 + (mem22/Em22)**2 + (mem22/Em22m)**2
    return  1 - 0.5*norm_c


def Vp00(delta, eta):
    return -(9*delta)/32

def Vp10(delta, eta):
    return 1/4

def Vm01(delta, eta):
    return -3/8

def Vm11(delta, eta): 
    return (3*delta)/64

def Vm21(delta, eta):
    return -(3*np.sqrt(2))/8

def Vp01(delta, eta): 
    return (3*eta*delta)/16

def Vp11(delta, eta):
    return -(3*eta)/2 

def Vm02(delta, eta):
    return - 9/64

def Vm12(delta, eta):
    return -(9*delta)/256

def Vm22(delta, eta):
    return -(9*np.sqrt(2))/64

Vp00_values = Vp00(delta, eta_list)
Vp10_values = Vp10(delta, eta_list)
Vm01_values = Vm01(delta, eta_list)
Vm11_values = Vm11(delta, eta_list)
Vm21_values = Vm21(delta, eta_list)
Vp01_values = Vp01(delta, eta_list)
Vp11_values = Vp11(delta, eta_list)
Vm02_values = Vm02(delta, eta_list)
Vm12_values = Vm12(delta, eta_list)
Vm22_values = Vm22(delta, eta_list)

Ep00_list = energy_0(1, 0, 0, eta_list)
Ep10_list = energy_0(1, 1, 0, eta_list)
Em01_list = energy_0(-1, 0, 1, eta_list)
Em01m_list = energy_0(-1, 0, -1, eta_list)
Em11_list = energy_0(-1, 1, 1, eta_list)
Em11m_list = energy_0(-1, 1, -1, eta_list)
Em21_list = energy_0(-1, 2, 1, eta_list)
Em21m_list = energy_0(-1, 2, -1, eta_list)
Ep01_list = energy_0(1, 0, 1, eta_list)
Ep01m_list = energy_0(1, 0, -1, eta_list)
Ep11_list = energy_0(1, 1, 1, eta_list)
Ep11m_list = energy_0(1, 1, -1, eta_list)
Em02_list = energy_0(-1, 0, 2, eta_list)
Em02m_list = energy_0(-1, 0, -2, eta_list)
Em12_list = energy_0(-1, 1, 2, eta_list)
Em12m_list = energy_0(-1, 1, -2, eta_list)
Em22_list = energy_0(-1, 2, 2, eta_list)
Em22m_list = energy_0(-1, 2, -2, eta_list)

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['font.size'] = 24


#fig = plt.figure(figsize=(7, 5))
#ax = fig.add_subplot(111)

energy_2_vardx = energy_2(delta_list, eta)
energy_2_vareta = energy_2(delta, eta_list)

#ax.set_xlim(3.5, 11.5)
#ax.set_xlim(0, 100)
#ax.set_ylim(-30010, -29990)
#ax.set_ylim(0, 1)

#ax.plot(delta_list, energy_2list)

# ax.plot(eta_list, Em01_list)
# ax.plot(eta_list, Em01m_list)
# ax.plot(eta_list, Em11_list)
# ax.plot(eta_list, Em11m_list)
# ax.plot(eta_list, Em21_list)
# ax.plot(eta_list, Em21m_list)
# ax.plot(eta_list, Em02_list)
# ax.plot(eta_list, Em02m_list)
# ax.plot(eta_list, Em12_list)
# ax.plot(eta_list, Em12m_list)
# ax.plot(eta_list, Em22_list)
# ax.plot(eta_list, Em22m_list)

#def energy_test(delta, eta):
#    return (((3*eta*delta)/16)**2)*((4*E)/(4*E**2 - eta**2))



#norm_m00_dx = norm_m00(delta_list, eta)
#norm_m00_eta = norm_m00(delta, eta_list)

#ax.plot(delta_list, energy_2_list)

#plt.tight_layout()
#plt.show()

