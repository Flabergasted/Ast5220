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

load = np.loadtxt("cells.txt")
ell     = load[:,0]
C_ell   = load[:,1]

print("C_ell")
fig, ax = plt.subplots()
fig.set_size_inches(10, 7.5)
ax.set_title(r"$\frac{\ell(\ell+1)C_\ell}{2\pi}$")
ax.plot(ell, (ell*(ell+1)*C_ell)/(2*np.pi), label = "with Peebles")
#ax.plot(x[Xe_Saha_of_x>1e-5], Xe_Saha_of_x[Xe_Saha_of_x>1e-5], ls = '--', label = "Saha only")
#ax.semilogy()
ax.loglog()
#ax.set_xlim(xmin = -12, xmax = 0)
#ax.set_ylim(ymin = 1e-4, ymax = 10)
ax.legend()
fig.savefig("C_ell")