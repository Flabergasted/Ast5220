import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from astropy import constants

Mpc = 1000 * constants.kpc.value
yr = const.year
kyr = 1e3*yr
Gyr = 1e9*yr
OmegaK0 = 0
OmegaB0 = 0.05
OmegaCDM0 = 0.267
Neff = 3.046
TCMB0 = 2.7255
H0 = 0.67*100
OmegaR0 = 2 * np.pi**2/30 * (const.k*TCMB0)**4/(const.hbar**3 * const.c**5) * (8*np.pi*const.G)/(3*H0)
OmegaNu0 = Neff * 7/8 * (4/11)**(4/3) * OmegaR0
OmegaLambda0 = 1 - (OmegaK0 + OmegaB0 + OmegaCDM0 + OmegaR0 + OmegaNu0)
Omega_rel0 = OmegaR0 + OmegaNu0
Omega_mass0 = OmegaB0 + OmegaCDM0

load = np.loadtxt("perturbations_k0.1.txt")
x           = load[:,0]
Theta0_1    = load[:,1]
Theta1_1    = load[:,2]
Theta2_1    = load[:,3]
Phi_1       = load[:,4]
Psi_1       = load[:,5]
Pi_1        = load[:,6]
delta_cdm_1 = load[:,7]
delta_b_1   = load[:,8]
v_cdm_1     = load[:,9]
v_b_1       = load[:,10]

load = np.loadtxt("perturbations_k0.01.txt")
#x_01        = load[:,0]
Theta0_01       = load[:,1]
Theta1_01       = load[:,2]
Theta2_01       = load[:,3]
Phi_01          = load[:,4]
Psi_01          = load[:,5]
Pi_01           = load[:,6]
delta_cdm_01    = load[:,7]
delta_b_01      = load[:,8]
v_cdm_01        = load[:,9]
v_b_01          = load[:,10]

load = np.loadtxt("perturbations_k0.001.txt")
#x_001        = load[:,0]
Theta0_001      = load[:,1]
Theta1_001      = load[:,2]
Theta2_001      = load[:,3]
Phi_001         = load[:,4]
Psi_001         = load[:,5]
Pi_001          = load[:,6]
delta_cdm_001   = load[:,7]
delta_b_001     = load[:,8]
v_cdm_001       = load[:,9]
v_b_001         = load[:,10]

plt.rcParams['font.size'] = 20

print("Density perturbations")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_title(r"Density perturbations")
ax.plot(x, np.abs(delta_cdm_1), color = "C2", label = r"k = 0.1 $Mpc^{-1}$")
ax.plot(x, np.abs(delta_cdm_01), color = "C1", label = r"k = 0.01 $Mpc^{-1}$")
ax.plot(x, np.abs(delta_cdm_001), color = "C0", label = r"k = 0.001 $Mpc^{-1}$")
ax.plot(x, np.abs(delta_b_1), color = "C2", ls = "--")
ax.plot(x, np.abs(delta_b_01), color = "C1", ls = "--")
ax.plot(x, np.abs(delta_b_001), color = "C0", ls = "--")
#"""
ax.plot(x, 4*Theta0_1, color = "C2", ls = "dotted")
ax.plot(x, 4*Theta0_01, color = "C1", ls = "dotted")
ax.plot(x, 4*Theta0_001, color = "C0", ls = "dotted")
#"""
ax.set_ylim(ymin = 1e-1, ymax = 1e5)
ax.set_xlabel("x")
ax.set_ylabel("$\delta_b$, $\delta_{CDM}$ and $\delta_{\gamma}$")
ax.semilogy()
ax.legend()
fig.savefig("density_perturbs")


print("Velocity perturbations")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_title(r"Velocity perturbations")
print("cdm full line")
ax.plot(x, np.abs(v_cdm_1), color = "C2", label = r"k = 0.1 $Mpc^{-1}$")
ax.plot(x, np.abs(v_cdm_01), color = "C1", label = r"k = 0.01 $Mpc^{-1}$")
ax.plot(x, np.abs(v_cdm_001), color = "C0", label = r"k = 0.001 $Mpc^{-1}$")
print("baryons --")
ax.plot(x, np.abs(v_b_1), color = "C2", ls = "--")
ax.plot(x, np.abs(v_b_01), color = "C1", ls = "--")
ax.plot(x, np.abs(v_b_001), color = "C0", ls = "--")
#"""
print("photons ..")
ax.plot(x, -3*Theta1_1, color = "C2", ls = "-.")
ax.plot(x, -3*Theta1_01, color = "C1", ls = "-.")
ax.plot(x, -3*Theta1_001, color = "C0", ls = "-.")
#"""
ax.semilogy()
ax.set_xlabel("x")
ax.set_ylabel(r"$v_b$, $v_{CDM}$, $v_\gamma$")
ax.legend()
fig.savefig("velocity_perturbs")


print("Gravitational Potential")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_title(r"Gravitational Potential")
ax.plot(x, Phi_1, label = r"k = 0.1 $Mpc^{-1}$")
ax.plot(x, Phi_01, label = r"k = 0.01 $Mpc^{-1}$")
ax.plot(x, Phi_001, label = r"k = 0.001 $Mpc^{-1}$")
ax.set_xlabel("x")
ax.set_ylabel(r"$\Phi$")
ax.legend()
fig.savefig("grav_potential")

print("Gravitational Potential + Psi")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_title(r"Gravitational Potential + Curvature")
ax.plot(x, Phi_1 + Psi_1, color = "C0", label = r"k = 0.1 $Mpc^{-1}$")
ax.plot(x, Phi_01 + Psi_01, color = "C1", label = r"k = 0.01 $Mpc^{-1}$")
ax.plot(x, Phi_001 + Psi_001, color = "C2", label = r"k = 0.001 $Mpc^{-1}$")
#ax.axvline(-8.3, color = "k", ls="--")
ax.set_xlabel("x")
ax.set_ylabel(r"$\Phi + \Psi$")
ax.legend()
fig.show()
fig.savefig("grav_potential_curvature")


print("Monopole")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_title(r"Monopole")
ax.plot(x, Theta0_001, label = r"k = 0.001 $Mpc^{-1}$")
ax.plot(x, Theta0_01, label = r"k = 0.01 $Mpc^{-1}$")
ax.plot(x, Theta0_1, label = r"k = 0.1 $Mpc^{-1}$")
ax.set_xlabel("x")
ax.set_ylabel(r"$\Theta_0$")
ax.legend()
fig.savefig("monopole")

print("Dipole")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_title("Dipole")
ax.plot(x, Theta1_001, label = r"k = 0.001 $Mpc^{-1}$")
ax.plot(x, Theta1_01, label = r"k = 0.01 $Mpc^{-1}$")
ax.plot(x, Theta1_1, label = r"k = 0.1 $Mpc^{-1}$")
ax.set_xlabel("x")
ax.set_ylabel(r"$\Theta_1$")
ax.legend()
fig.savefig("dipole")

plt.rcParams['font.size'] = 15

print("Quadrapole")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_title("Quadrapole")
ax.plot(x, Theta2_001, label = r"k = 0.001 $Mpc^{-1}$")
ax.plot(x, Theta2_01, label = r"k = 0.01 $Mpc^{-1}$")
ax.plot(x, Theta2_1, label = r"k = 0.1 $Mpc^{-1}$")
ax.set_xlabel("x")
ax.set_ylabel(r"$\Theta_2$")
ax.legend()
fig.savefig("quadrapole")