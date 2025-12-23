import numpy as np
import matplotlib.pyplot as plt
import os

# -------------------
# Parameters & ranges
# -------------------
E = np.logspace(2, 8, 500)  # Energy range: 10^2 to 10^8 GeV
E0 = 5e4                    # Reference energy [GeV]
phi0_values = [2e-14, 2e-13, 2e-12]
gamma_values = [1, 1.5, 2, 2.5, 3]
Emin, Emax = 5e4, 5e5        # Bounds for the bounded flux

E = np.logspace(0, 8, 500)  # Energy range: 10^2 to 10^8 GeV
E0 = 1                    # Reference energy [GeV]
phi0_values = [2e-6, 2e-5, 2e-3]
gamma_values = [1, 1.5, 2, 2.5, 3]
Emin, Emax = 5e4, 5e5        # Bounds for the bounded flux

# Output folder
out_dir = "flux_plots"
os.makedirs(out_dir, exist_ok=True)

# -------------------
# Flux definitions
# -------------------
def phi_classic(E, phi0, gamma):
    """Classic unbounded power-law flux."""
    return phi0 * (E / E0)**(-gamma)

def phi_bounded(E, phi0, gamma):
    """Bounded power-law flux: zero outside [Emin, Emax]."""
    mask = (E >= Emin) & (E <= Emax)
    flux = np.zeros_like(E)
    flux[mask] = phi_classic(E[mask], phi0, gamma)
    return flux

# -------------------
# Generate and save each plot
# -------------------
for gamma in gamma_values:
    for phi0 in phi0_values:
        flux_class = phi_classic(E, phi0, gamma)
        #flux_bound = phi_bounded(E, phi0, gamma)

        fig, ax = plt.subplots(figsize=(5, 4))
        ax.fill_between(E, flux_class, alpha=1, color='darkblue' , label=r'$\Phi_{\mathrm{classic}}$')
        #ax.fill_between(E, flux_bound, alpha=1, color='darkorange', label=r'$\Phi_{\mathrm{bounded}}$')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(1e2, 1e8)
        ax.set_ylim(1e-17, 1e-10)
        #ax.set_title(fr'$\gamma={gamma}, \Phi_0={phi0:.0e}$')
        ax.set_xlabel('Energy [GeV]')
        ax.set_ylabel(r'Flux $\Phi(E)$ [GeV$^{-1}$ m$^{-2}$ s$^{-1}$]')
        ax.grid(True, ls="--", alpha=0.5)
        ax.legend()

        # Filenames
        filename_base = f"unbounded_flux_gamma{gamma}_phi{phi0:.0e}"
        png_path = os.path.join(out_dir, f"{filename_base}.png")
        pdf_path = os.path.join(out_dir, f"{filename_base}.pdf")

        # Save files
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        plt.close(fig)

print(f"Saved plots in '{out_dir}' folder.")

import numpy as np
import matplotlib.pyplot as plt
import os

# -------------------
# Parameters & ranges
# -------------------
  # Bounds for the bounded flux

# Output folder
out_dir = "flux_plots"
os.makedirs(out_dir, exist_ok=True)

# -------------------
# Flux definitions
# -------------------
def phi_classic(E, phi0, gamma):
    """Classic unbounded power-law flux."""
    return phi0 * (E / E0)**(-gamma)

def phi_bounded(E, phi0, gamma):
    """Bounded power-law flux: zero outside [Emin, Emax]."""
    mask = (E >= Emin) & (E <= Emax)
    flux = np.zeros_like(E)
    flux[mask] = phi_classic(E[mask], phi0, gamma)
    return flux

# -------------------
# Generate and save plots for each phi0
# -------------------
for phi0 in phi0_values:
    fig, ax = plt.subplots(figsize=(6, 5))

    # Example list of colors (can be any valid Matplotlib colors)
    colors = ['teal','steelblue','darkblue', 'purple', 'firebrick'] 
    alphas = [1,0.9,0.8,0.7,0.6]

    # Sort gammas to plot from lowest to highest
    gammas_sorted = sorted(gamma_values)

    for i, gamma in enumerate(gammas_sorted):
        #flux_bound = phi_bounded(E, phi0, gamma)
        flux_bound = phi_classic(E, phi0, gamma)
        # Use color from the list (wrap around if more gammas than colors)
        color = colors[i % len(colors)]
        alpha = alphas[i % len(alphas)]
        
        ax.plot(E, flux_bound, lw=2, label=fr'$\gamma={gamma}$', color=color)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e0, 1e8)
    ax.set_ylim(1e-16, 1e-2)
    ax.set_xlabel('Energy [GeV]')
    ax.set_ylabel(r'Flux $\Phi(E)$ [GeV$^{-1}$ m$^{-2}$ s$^{-1}$]')
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.legend()
    ax.set_title(fr'Unbounded Pwer Law Flux for $\Phi_0={phi0:.0e}$')

    # Save files
    filename_base = f"flux_unbounded_phi{phi0:.0e}"
    png_path = os.path.join(out_dir, f"{filename_base}.png")
    pdf_path = os.path.join(out_dir, f"{filename_base}.pdf")
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path)
    plt.close(fig)

print(f"Saved plots in '{out_dir}' folder.")

