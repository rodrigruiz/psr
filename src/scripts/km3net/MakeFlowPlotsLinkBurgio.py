import numpy as np
import matplotlib.pyplot as plt
import os

# -------------------
# Output folder
# -------------------
out_dir = "flux_plots"
os.makedirs(out_dir, exist_ok=True)

# -------------------
# Flux definitions
# -------------------
def phi_classic(E, phi0, gamma, E0):
    """Classic unbounded power-law flux."""
    return phi0 * (E / E0)**(-gamma)

def phi_bounded(E, phi0, gamma, E0, Emin, Emax):
    """Bounded power-law flux: zero outside [Emin, Emax]."""
    mask = (E >= Emin) & (E <= Emax)
    flux = np.zeros_like(E)
    flux[mask] = phi_classic(E[mask], phi0, gamma, E0)
    return flux

    import numpy as np

def phi_softbounded(E, phi0, gamma, E0, Emin, Emax, delta = 50, factor = 2, xoffset = 0):
    """
    Power-law flux with soft lower cutoff.
    
    E     : energy array
    phi0  : normalization at E0
    gamma : spectral index
    E0    : reference energy
    Emin  : characteristic lower cutoff energy
    Emax  : sharp upper cutoff
    delta : softening scale of the cutoff (controls smoothness)
    """
    # Core power-law
    flux = factor* phi0 * ((E+xoffset) / (E0-delta))**(-gamma)
    
    # Soft turn-on near Emin
    soft = 1 - np.exp(-((E+xoffset) - (Emin-delta)) / delta)
    
    # Upper cutoff
    mask =(E+xoffset) <= Emax
    return flux * soft * mask


import numpy as np

def phi_softlower(E, phi0, gamma, E0, Emin, Emax, delta, eta=0.0):
    """
    Power-law with soft lower roll-off anchored at E0 (with E0 == Emin),
    and a (still hard) upper cutoff at Emax.

    Parameters
    ----------
    E : array_like
        Energies.
    phi0 : float
        Normalization: phi(E0) = phi0.
    gamma : float
        Spectral index of the power law.
    E0 : float
        Reference energy; also the soft-cut anchor (Emin). Peak is at (E0, phi0).
    Emax : float
        Hard upper cutoff energy (set flux=0 above this).
    delta : float
        Softening scale (in energy units) controlling how quickly the flux
        rolls off to the left of E0. Larger = gentler roll-off.
    eta : float, optional (default 0.0)
        Extra left-side power suppression. If gamma>0 and you want the low-E
        tail to vanish (and avoid any divergence as E→0), set eta > gamma.

    Notes
    -----
    - The left soft factor is:
        soft_left(E) = exp(- max(E0 - E, 0) / delta) * [(E/E0)^eta for E<E0, else 1]
      So:
        * E >= E0  → soft_left = 1  (so phi(E0)=phi0, untouched to the right)
        * E <  E0  → exponential falloff over scale `delta`, optionally steeper by (E/E0)^eta
    """
    E = np.asarray(E)

    # Base power-law (assumes E>0 on your grid)
    flux = phi0 * (E / E0)**(-gamma)

    # Exponential roll-off only for E<E0 (equals 1 at and above E0)
    soft_exp = np.exp(-np.clip(E0 - E, 0.0, None) / delta)

    # Optional extra power suppression on the left to tame low-E behavior
    if eta > 0.0:
        soft_pow = np.where(E < E0, (E / E0)**eta, 1.0)
        soft = soft_exp * soft_pow
    else:
        soft = soft_exp

    # Apply (still hard) upper cutoff
    flux = flux * soft
    flux = np.where(E <= Emax, flux, 0.0)

    return flux


# -------------------
# Parameters set 1
# -------------------
params1 = {
    "E": np.logspace(4, 8, 500),
    "E0": 1.95e5,
    "phi0": 8e-15,
    "gamma": 2.1,
    "Emin": 2e5,
    "Emax": 2e6,
    'delta': 30000,
    'xoffset':15000 
}

# -------------------
# Parameters set 2
# -------------------
params2 = {
    "E": np.logspace(4, 8, 500),
    "E0": 1.3e5,
    "phi0": 5e-15,
    "gamma": 2.6,
    "Emin": 1.3e5,
    "Emax": 2e6,
    'delta': 15000,
    'xoffset':0  
}

# -------------------
# Plot both in same figure
# -------------------
fig, ax = plt.subplots(figsize=(6, 5))

# Style options
colors = ["k", "k"]
labels = [fr"Linear origin:        $\gamma={params1['gamma']}$   ,  $\Phi_0={params1['phi0']:.0e}$ ,  $E_{0}={params1['E0']:.0e}$ ,  $E_{{max}}={params1['Emax']:.0e}$",
          fr"Quadratic origin:  $\gamma={params2['gamma']}$ ,  $\Phi_0={params2['phi0']:.0e}$ ,  $E_{0}={params2['E0']:.0e}$ ,  $E_{{max}}={params2['Emax']:.0e}$"]

# Compute fluxes (choose classic or bounded)
#flux1 = phi_bounded(params1["E"], params1["phi0"], params1["gamma"], params1["E0"])
#flux2 = phi_bounded(params2["E"], params2["phi0"], params2["gamma"], params2["E0"])
#flux1 = phi_bounded(**params1)
#flux2 = phi_bounded(**params2)
flux1 = phi_softbounded(**params1)
flux2 = phi_softbounded(**params2)
#flux1 = phi_softlower(**params1)
#flux2 = phi_softlower(**params2)

ax.plot(params1["E"], flux1, lw=2, color=colors[0], label=labels[0])
ax.plot(params2["E"], flux2, lw=2, ls='--', color=colors[1], label=labels[1])

# Axes scaling and labels
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e4, 1e7)
ax.set_ylim(1e-16, 1e-12)
ax.set_xlabel("Energy [GeV]")
ax.set_ylabel(r"Flux $\Phi(E)$ [GeV$^{-1}$ m$^{-2}$ s$^{-1}$]")
ax.grid(True, which="both", ls="--", alpha=0.5)
ax.legend(fontsize=8,loc="upper center")
ax.set_title("Remodeled Link&Burgio Flux Prediction for Vela")

# Save
filename_base = "flux_vela_linkburgionew"
plt.savefig(os.path.join(out_dir, f"{filename_base}.png"), dpi=300)
plt.savefig(os.path.join(out_dir, f"{filename_base}.pdf"))
plt.close(fig)

print(f"Saved plots in '{out_dir}' folder.")
