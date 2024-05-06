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

load = np.loadtxt("recombination.txt")
x                   = load[:,0]
Xe_of_x             = load[:,1]
ne_of_x             = load[:,2]
tau_of_x            = load[:,3]
dtaudx_of_x         = load[:,4]
ddtauddx_of_x       = load[:,5]
g_tilde_of_x        = load[:,6]
dgdx_tilde_of_x     = load[:,7]
ddgddx_tilde_of_x   = load[:,8]
t_of_x              = load[:,9]
s_of_x              = load[:,10]
Xe_Saha_of_x        = load[:,11]

plt.rcParams['font.size'] = 15

print("Xe_of_x")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_title(r"$X_e (x)$")
ax.plot(x, Xe_of_x, label = "with Peebles")
print(Xe_Saha_of_x)
ax.plot(x[Xe_Saha_of_x>1e-5], Xe_Saha_of_x[Xe_Saha_of_x>1e-5], ls = '--', label = "Saha only")
ax.semilogy()
ax.set_xlim(xmin = -12, xmax = 0)
ax.set_ylim(ymin = 1e-4, ymax = 10)
ax.legend()
fig.savefig("Xe_of_x")

print("tau_of_x")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_title(r"$\tau (x)$")
ax.plot(x, tau_of_x, label=r"$\tau(x)$")
ax.plot(x, -dtaudx_of_x, label=r"$-\tau'(x)$")
ax.plot(x, ddtauddx_of_x, label=r"$\tau''(x)$")
ax.legend()
ax.semilogy()
ax.set_xlim(xmin = -12, xmax = 0)
#ax.set_ylim(ymin = 1e-8, ymax = 1e8)
fig.savefig("tau_of_x")

print("g_tilde_of_x")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_title(r"$\tilde{g} (x)$")
ax.plot(x, g_tilde_of_x/np.max(g_tilde_of_x), label=r"$\tilde{g}(x)$")
ax.plot(x, dgdx_tilde_of_x/np.max(dgdx_tilde_of_x), ls="-.", label=r"$\tilde{g}'(x)$")
ax.plot(x, ddgddx_tilde_of_x/np.max(ddgddx_tilde_of_x), ls="--", label=r"$\tilde{g}''(x)$")
ax.set_xlim(xmin = -8, xmax = -6)
ax.legend()
#ax.set_ylim(ymin = 1e-8, ymax = 1e8)
fig.savefig("g_tilde_of_x")

x_Xe = x[np.argmin(np.abs(Xe_of_x-0.1))]
x_Saha_Xe = x[np.argmin(np.abs(Xe_Saha_of_x-0.1))]
x_tau = x[np.argmin(np.abs(tau_of_x-1))]
x_g_tilde = x[np.argmax(g_tilde_of_x)]

print("Time of recombination & last scattering")
print("")
print("For x:")
print("Recombination, Xe = 0.1:", x_Xe)
print("Recombination, Saha:", x_Saha_Xe)
print("Last scattering, tau(x) = 1:", x_tau)
print("Last scattering, max(g_tilde(x)):",x_g_tilde)
print("")
print("For a:")
print("Recombination, Xe = 0.1:", np.exp(x_Xe))
print("Recombination, Saha:", np.exp(x_Saha_Xe))
print("Last scattering, tau(x) = 1:", np.exp(x_tau))
print("Last scattering, max(g_tilde(x)):", np.exp(x_g_tilde))
print("")
print("For z:")
print("Recombination, Xe = 0.1:", 1/np.exp(x_Xe) - 1)
print("Recombination, Saha:", 1/np.exp(x_Saha_Xe) - 1)
print("Last scattering, tau(x) = 1:", 1/np.exp(x_tau) - 1)
print("Last scattering, max(g_tilde(x)):", 1/np.exp(x_g_tilde) - 1)
print("")
print("For t:")
print("Recombination, Xe = 0.1:", t_of_x[np.argmin(np.abs(Xe_of_x-0.1))]/kyr," kyr")
print("Recombination, Saha:", t_of_x[np.argmin(np.abs(Xe_Saha_of_x-0.1))]/kyr," kyr")
print("Last scattering, tau(x) = 1:", t_of_x[np.argmin(np.abs(tau_of_x-1))]/kyr," kyr")
print("Last scattering, max(g_tilde(x)):",t_of_x[np.argmax(g_tilde_of_x)]/kyr," kyr")
print("")
print("Free electrons today", Xe_of_x[-1])
print("Sound horizon: ",s_of_x[np.argmin(np.abs(x-x_tau))]/Mpc, " Mpc")
