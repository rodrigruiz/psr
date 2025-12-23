#!/usr/bin/env python3
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
import os
import csv

# ---- CREATE OUTPUT FOLDER ----
output_folder = "output"
os.makedirs(output_folder, exist_ok=True)

# ---- PARAMETERS ----
R = 1e6           # stellar radius in meters (10 km)
L = 1.0           # acceleration length in units of R
xw = 0.1          # resonance width
sigma0 = 5e-28    # cm^2 -> m^2
sigma0 *= 1e-4    # convert to m^2
phi0 = 1.0        # normalization
z_min, z_max = 0, 5
z_vals = np.linspace(z_min, z_max, 200)
x_vals = np.linspace(1.1, 50, 200)
E_nu = 0.05 * x_vals  # neutrino energy fraction

# ---- FUNCTIONS ----
def theta_h(z):
    return np.arcsin(1 / (1 + z))

def theta_0(x):
    cos_theta0 = 1 - 1/(x + xw)
    return np.arccos(np.clip(cos_theta0, -1, 1))

def dn_dtheta(theta, z):
    alpha = np.arcsin((1 + z) * np.sin(theta)) - theta
    d_alpha_dtheta = (1 + z) * np.cos(theta) / np.sqrt(1 - (1 + z)**2 * np.sin(theta)**2) - 1
    return np.sin(theta)**2 * np.sin(alpha) * np.cos(alpha + theta) * d_alpha_dtheta

def sigma(x, theta):
    xT = 1 / (1 - np.cos(theta))
    return sigma0 * np.exp(-(x - xT)**2 / (2 * xw**2)) if x > xT - xw else 0.0

def integrand(theta, x, z):
    return dn_dtheta(theta, z) * sigma(x, theta)

def dPc_dx(x, gamma):
    total = 0
    for z in z_vals:
        theta0 = theta_0(x)
        thetah = theta_h(z)
        if theta0 >= thetah:
            continue
        result, _ = quad(integrand, theta0, thetah, args=(x, z))
        dz_dx = L / gamma * x**(1/gamma - 1)
        total += result * dz_dx
    return R * total

# ---- COMPUTE FLUX ----
fluxes = {}
for gamma in [1, 2]:
    flux = phi0 * np.array([dPc_dx(x, gamma) for x in x_vals])
    fluxes[gamma] = flux

# ---- PLOT ----
plt.figure(figsize=(8,5))
plt.loglog(E_nu, fluxes[1], label="Linear acceleration (γ=1)")
plt.loglog(E_nu, fluxes[2], label="Quadratic acceleration (γ=2)")
plt.xlabel(r'$E_\nu$ [arbitrary units]')
plt.ylabel(r'$d\phi_\nu/dE_\nu$')
plt.title('Neutrino flux from pulsar (full z-dependence)')
plt.grid(True, which="both", ls="--", lw=0.5)
plt.legend()
plot_file = os.path.join(output_folder, "neutrino_flux.png")
plt.savefig(plot_file)
print(f"Plot saved to {plot_file}")
plt.close()

# ---- SAVE FLUX TO CSV ----
csv_file = os.path.join(output_folder, "neutrino_flux.csv")
with open(csv_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["E_nu"] + [f"flux_gamma{g}" for g in [1, 2]])
    for i in range(len(E_nu)):
        writer.writerow([E_nu[i], fluxes[1][i], fluxes[2][i]])
print(f"Flux values saved to {csv_file}")
