import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from astropy import constants

Mpc = 1000 * constants.kpc.value
Gyr = 10**9*365*24*60*60

load = np.loadtxt("cosmology.txt")
x = load[:,0]
eta = load[:,1]
t_of_x = load[:,2]
Hp_of_x = load[:,3]
dHpdx = load[:,4]
ddHpddx = load[:,5]
OmegaB_of_x = load[:,6]
OmegaCDM_of_x = load[:,7]
OmegaLambda_of_x = load[:,8]
OmegaR_of_x = load[:,9]
OmegaNu_of_x = load[:,10]
OmegaK_of_x = load[:,11]
luminosity_distance = load[:,12]
comoving_distance = load[:,13]

print(f"the age of the universe x = {x[-1]}, eta (Gyr) = {eta[-1]/Gyr}, t (Gyr) = {t_of_x[-1]/Gyr}")

supernova = np.loadtxt("results_supernovafitting.txt")
chi = supernova[:,0]
h = supernova[:,1]
OmegaM = supernova[:,2]
OmegaK = supernova[:,3]
#2 sigma = 8.02
sigma1 = chi[(chi - np.min(chi) < 3.53)]
passed_args = np.argwhere((np.abs(chi - np.min(chi)) < 3.53))
passed_args2 = np.argwhere((np.abs(chi - np.min(chi)) < 8.02))

OmegaM_sigma = OmegaM[passed_args]
OmegaK_sigma = OmegaK[passed_args]
OmegaLambda_sigma = 1 - OmegaK_sigma - OmegaM_sigma
OmegaM_sigma2 = OmegaM[passed_args2]
OmegaK_sigma2 = OmegaK[passed_args2]
OmegaLambda_sigma2 = 1 - OmegaK_sigma2 - OmegaM_sigma2
#print(chi[passed_args] - np.min(chi))
print("Begin plot creation")

print("1 sigma constraints from MCMC in the OmegaLambda OmegaM plane")
fig, ax = plt.subplots()
ax.set_title(r"$1\sigma$ deviation in the $\Omega_\Lambda$-$\Omega_m$-plane")
ax.scatter(OmegaM_sigma2, OmegaLambda_sigma2, label = r"$2\sigma$")
ax.scatter(OmegaM_sigma, OmegaLambda_sigma, label = r"$1\sigma$")
ax.plot(OmegaM, 1 - OmegaM, ls = "-.", color = "k", label = "flat universe")
ax.set_xlim(xmin = 0, xmax = 1)
ax.set_ylim(ymin = 0, ymax = 1.4)
ax.set_xlabel(r"$\Omega_M$")
ax.set_ylabel(r"$\Omega_\Lambda$")
ax.legend()
fig.savefig("sigma_M_Lambda")
#print(supernova[(supernova[:,0] - np.min(supernova[:,0]) < 3.53),0])

mass = OmegaB_of_x + OmegaCDM_of_x
relativistic = OmegaR_of_x + OmegaNu_of_x

print("Omegas")
fig, ax = plt.subplots()
ax.set_title(r"Evolution of $\Omega_i(x)$")
ax.plot(x, mass, label = r"$\Omega_m + \Omega_{CDM}$")
ax.plot(x, relativistic, label = r"$\Omega_\gamma + \Omega_\nu$")
ax.plot(x, OmegaLambda_of_x, label = r"$\Omega_{\Lambda}$")
ax.set_xlabel(r"$x = ln(a)$")
ax.set_ylabel("Mass fraction in universe")
ax.axvline(np.log(2.94*1e-4), ls = "--", color="k", label = "Matter-Radiation equality")
ax.axvline(np.log(0.775), ls = "--", color="r", label = "Matter-Dark energy equality")
ax.set_ylim(ymin = -0.05, ymax = 1.05)
ax.legend()
ax.grid()
fig.savefig("Omegas_M1.png")

fig, ax = plt.subplots()
ax.set_title(r"Deviation in $\sum\Omega_i$")
ax.plot(x, np.abs(1 - (mass+relativistic+OmegaLambda_of_x)))
ax.set_xlabel(r"$x = ln(a)$")
ax.set_ylabel("Deviation from 1")
ax.semilogy()
fig.savefig("Omega_deviation_from_1.png")

print("eta/c vs t")
fig, ax = plt.subplots()
ax.set_title(r"$\frac{\eta (x)}{c}$ vs $t$ in Gyr")
ax.plot(x, t_of_x/Gyr, label = "t")
ax.plot(x, eta/(const.c*Gyr), label = r"$\eta$")
ax.set_xlabel("x = log(a)")
ax.set_ylabel("Age in Gyrs")
ax.semilogy()
ax.legend()
ax.grid()
fig.savefig("eta_of_t.png")

print("eta(x)")
fig, ax = plt.subplots()
ax.set_title(r"$\eta (x)$ (Mpc)")
ax.plot(x, eta/Mpc)
ax.semilogy()
ax.set_xlabel("x = log(a)")
ax.set_ylabel(r"$\eta (x)$ in Mpc")
ax.grid()
fig.savefig("eta_of_x.png")

print("Hp")
fig, ax = plt.subplots()
ax.set_title(r"$\mathcal{H} (\frac{100 km/s}{Mpc})$")
ax.plot(x, Hp_of_x*Mpc/(100*1000))
ax.semilogy()
ax.set_xlabel("x = log(a)")
ax.set_ylabel(r"$\mathcal{H} (x)$")
ax.grid()
fig.savefig("Hp.png")

print("eta*Hp/c")
fig, ax = plt.subplots()
ax.set_title(r"$\frac{\eta(x) \mathcal{H}(x)}{c}$")
ax.plot(x, eta*Hp_of_x/const.c)
ax.set_xlabel("x = log(a)")
ax.set_ylabel(r"$\eta (x) \mathcal{H} (x) / c$")
ax.grid()
fig.savefig("eta_Hp_over_c.png")

fig, ax = plt.subplots()
print("1/Hp* dHp/dx")
ax.set_title(r"$\frac{1}{\mathcal{H}(x)} \frac{\mathrm{d}\mathcal{H}(x)}{\mathrm{d}x}$")
ax.plot(x, dHpdx/Hp_of_x)
ax.set_xlabel("x = log(a)")
ax.set_ylabel("$\mathcal{H}(x)^{-1} \mathcal{H}'(x)$")
#ax.hline()
ax.grid()
fig.savefig("dHpdx")

fig, ax = plt.subplots()
print("1/Hp* d^2Hp/dx^2")
ax.set_title(r"$\frac{1}{\mathcal{H}(x)} \frac{\mathrm{d}^2\mathcal{H}(x)}{\mathrm{d}x^2}$")
ax.plot(x, ddHpddx/Hp_of_x)
ax.set_xlabel("x = log(a)")
ax.set_ylabel(r"$\mathcal{H}(x)^{-1} \mathcal{H}''(x)$")
#ax.plot(x[90:], -1/2 * (0.67*100)**2 * (9 * 0.317)/((np.exp(x[90:]))**3), ls="--")
ax.grid()
fig.savefig("ddHpddx")

SN_data = np.loadtxt("./data/supernovadata.txt")
z = SN_data[:,0]
d_L_data = SN_data[:,1]
SN_errorbars = SN_data[:,2]

print("plot supernova data vs z here")
fig, ax = plt.subplots()
ax.set_title(r"$d_L$ supernova data vs fitting")
ax.errorbar(z, d_L_data, yerr=SN_errorbars, fmt=".", capsize=3, color="k", label = "data")
a_temp = np.exp(x)
z_temp = 1/a_temp - 1
ax.plot(z_temp, luminosity_distance/(1000*Mpc), label = "fitting")
ax.set_xlabel("z")
ax.set_ylabel(r"$d_L$ (Gpc)")
ax.set_xlim(xmin = 8*1e-3, xmax = 2)
ax.set_ylim(ymin = -0.5, ymax = 10)
ax.semilogx()
ax.grid()
ax.legend()
fig.savefig("supernova_data_plot.png")

def gaussian(input_param):
    return 1/(np.sqrt(2*np.pi)*np.std(input_param))*np.exp(-1/2 * (input_param-np.mean(input_param))**2/np.std(input_param)**2)

print("Poseterior PDF of H0")
fig, ax = plt.subplots(2,2, figsize = (10,10))
ax[0,0].set_title(r"Posterior for $H_0$")
ax[0,0].hist(h, bins = 100, density=True)
ax[0,0].set_xlim(xmax = 0.73)
ax[0,0].plot(np.sort(h),gaussian(np.sort(h, axis=None)))
ax[0,0].set_ylabel("Probability density")
#ax[0,0].vlines(0.67, ymin=0, ymax=1000, ls="--", color="k", label="Planck (2018) best-fit value")
#fig.savefig("H0_posterior_hist.png")

print("Poseterior PDF of OmegaM")
#fig, ax = plt.subplots()
ax[0,1].set_title(r"Posterior for $\Omega_M$")
ax[0,1].hist(OmegaM, bins = 40, density=True)
#ax[0,1].set_xlim(xmin = 0.66, xmax = 0.73)
ax[0,1].vlines(0.05 + 0.267, ymin=0, ymax=4, ls="--", color="k", label="Planck (2018) best-fit value")
ax[0,1].plot(np.sort(OmegaM),gaussian(np.sort(OmegaM, axis=None)))
ax[0,1].legend()
#fig.savefig("OmegaM_posterior_hist.png")

print("Poseterior PDF of OmegaK")
#fig, ax = plt.subplots()
ax[1,0].set_title(r"Posterior for $\Omega_K$")
ax[1,0].hist(OmegaK, bins = 40, density=True)
#ax[1,0].set_xlim(xmin = 0.66, xmax = 0.73)
ax[1,0].vlines(0, ymin=0, ymax=1.5, ls="--", color="k", label="Planck (2018) best-fit value")
ax[1,0].plot(np.sort(OmegaK),gaussian(np.sort(OmegaK, axis=None)))
ax[1,0].legend()
ax[1,0].set_xlabel("predicted values")
ax[1,0].set_ylabel("Probability density")
#fig.savefig("OmegaK_posterior_hist.png")

print("Poseterior PDF of OmegaLambda")
OmegaK0 = 0
OmegaB0 = 0.05
OmegaCDM0 = 0.267
Neff = 3.046
TCMB0 = 2.7255
H0 = 0.67*100
OmegaR0 = 2 * np.pi**2/30 * (const.k*TCMB0)**4/(const.hbar**3 * const.c**5) * (8*np.pi*const.G)/(3*H0)
OmegaNu0 = Neff * 7/8 * (4/11)**(4/3) * OmegaR0
OmegaLambda_best_fit = 1 - (OmegaK0 + OmegaB0 + OmegaCDM0 + OmegaR0 + OmegaNu0)

#fig, ax = plt.subplots()
ax[1,1].set_title(r"Posterior for $\Omega_\Lambda$")
ax[1,1].hist(1 - (OmegaK + OmegaM), bins = 40, density=True)
ax[1,1].vlines(OmegaLambda_best_fit, ymin=0, ymax=10, ls="--", color="k", label="Planck (2018) best-fit value")
ax[1,1].plot(np.sort(1 - (OmegaK + OmegaM)),gaussian(np.sort(1 - (OmegaK + OmegaM), axis=None)))
ax[1,1].set_ylim(ymax = 3)
ax[1,1].legend()
ax[1,1].set_xlabel("predicted values")
#fig.savefig("OmegaLambda_posterior_hist.png")
fig.savefig("Posterior_histograms.png")

print("complete")