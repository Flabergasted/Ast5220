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
h = 0.67
H0 = h*100
OmegaR0 = 2 * np.pi**2/30 * (const.k*TCMB0)**4/(const.hbar**3 * const.c**5) * (8*np.pi*const.G)/(3*H0)
OmegaNu0 = Neff * 7/8 * (4/11)**(4/3) * OmegaR0
OmegaLambda0 = 1 - (OmegaK0 + OmegaB0 + OmegaCDM0 + OmegaR0 + OmegaNu0)
Omega_rel0 = OmegaR0 + OmegaNu0
Omega_mass0 = OmegaB0 + OmegaCDM0

cosmo = np.loadtxt("cosmology.txt")
x                   = cosmo[:,0]
eta0                = cosmo[-1,1]
Hp_of_x             = cosmo[:,3]
MatterRadiationEq   = np.log(2.94*1e-4)
x_eq = np.argmin(np.abs(x-MatterRadiationEq))
k_eq = Hp_of_x[x_eq]/(const.c)

load = np.loadtxt("cells.txt")
ell         = load[:,0]
C_ell       = load[:,1]

data = np.loadtxt("planck_cell_low.txt")
ELL = data[:,0]
C_ELL = data[:,1]
ERR_UP = data[:,2]
ERR_DOWN = data[:,3]

errors = (ERR_DOWN, ERR_UP)

load2 = np.loadtxt("Pk_ell.txt")
k = load2[:,0]
P_m = load2[:,1]
theta = load2[:,2]
theta5 = load2[:,3]
theta50 = load2[:,4]
theta100 = load2[:,5]
theta500 =  load2[:,6]

plt.rcParams['font.size'] = 20

print("CMB power spectrum, Cell")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_title(r"The CMB Power spectrum $C_\ell$")
ax.set_ylabel(r"$\ell(\ell+1)C_\ell/2\pi$ $(\mu K)$")
ax.set_xlabel(r"Multipole $\ell$")
ax.plot(ell, C_ell, label = r"$C_\ell$")
ax.errorbar(ELL, C_ELL, errors, capsize = 3, marker = ".", ls = "none", color = "k")
ax.legend()
ax.set_ylim(ymin = 1e2, ymax=1e4)
plt.loglog()
fig.savefig("C_ell")

plt.rcParams['font.size'] = 15

print("Theta_ell_0")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_ylabel(r"$\Theta_0(k)^2$")
ax.plot(k*eta0, theta)
plt.semilogx()
ax.set_xlabel(r"$k\eta_0\approx\ell$")
fig.savefig("Theta_ell")

print("Theta_ell_5")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_ylabel(r"$\Theta_5(k)^2$")
ax.plot(k*eta0, theta5)
ax.axvline(5,ls='--')
plt.semilogx()
fig.savefig("Theta_ell_5")

print("Theta_ell_50")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_ylabel(r"$\Theta_{50}(k)^2$")
ax.plot(k*eta0, theta50)
ax.axvline(50,ls='--')
plt.semilogx()
ax.set_xlabel(r"$k\eta_0\approx\ell$")
fig.savefig("Theta_ell_50")

print("Theta_ell_100")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_ylabel(r"$\Theta_{100}(k)^2$")
ax.plot(k*eta0, theta100)
ax.axvline(100,ls='--')
plt.semilogx()
ax.set_xlabel(r"$k\eta_0\approx\ell$")
fig.savefig("Theta_ell_100")

print("Theta_ell_500")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_ylabel(r"$\Theta_{500}(k)^2$")
ax.plot(k*eta0, theta500)
ax.axvline(500,ls='--')
plt.semilogx()
ax.set_xlabel(r"$k\eta_0\approx\ell$")
fig.savefig("Theta_ell_500")

plt.rcParams['font.size'] = 20

print("Theta_ells_sq")
fig = plt.figure()
gs = fig.add_gridspec(3, hspace=0)
axs = gs.subplots(sharex=True, sharey=False)
fig.set_size_inches(10, 7.5)
fig.suptitle(r"Squared $\Theta_\ell$s")
fig.supxlabel(r"$k\eta_0\approx\ell$")
fig.supylabel(r"$|\Theta_\ell(k)|^2/k/(H_0^{-1})$")

axs[0].plot(k*eta0, np.abs(theta5)**2/(k*Mpc)/(H0**(-1)), label=r"$\Theta_5$")
axs[1].plot(k*eta0, np.abs(theta50)**2/(k*Mpc)/(H0**(-1)), label=r"$\Theta_{50}$")
axs[2].plot(k*eta0, np.abs(theta100)**2/(k*Mpc)/(H0**(-1)), label=r"$\Theta_{100}$")

axs[0].legend()
axs[1].legend()
axs[2].legend()
plt.xlim(xmin=0, xmax = 500)
plt.tight_layout()
#axs[0].set_xlim(xmin=0, xmax = 500)
#axs[1].set_xlim(xmin=0, xmax = 500)
#axs[2].set_xlim(xmin=0, xmax = 500)
fig.savefig("Theta_ells_sq")

print("Matter power spectrum today with the equality scale k_eq marked")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_title("Matter Power spectrum")
ax.set_ylabel(r"$P(k)$ $(Mpc/h)^3$")
ax.plot(k*(Mpc/h), P_m)
ax.axvline(k_eq*(Mpc/h),ls="--")
plt.loglog()
plt.grid()
ax.set_xlabel(r"$k (Mpc/h)$")
fig.savefig("matter_power_spec")

print(C_ell[np.argmax(C_ell)])