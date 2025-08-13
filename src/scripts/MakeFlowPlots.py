import numpy as np
import matplotlib.pyplot as plt
import os

# -------------------
# Parameters & ranges
# -------------------
E = np.logspace(2, 8, 50)  # Energy range: 10^2 to 10^8 GeV
E0 = 1e5                    # Reference energy [GeV]
phi0_values = [1e-11, 1e-10]
gamma_values = [1.5, 2.5]
Emin, Emax = 5e4, 5e6        # Bounds for the bounded flux

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
        flux_bound = phi_bounded(E, phi0, gamma)

        fig, ax = plt.subplots(figsize=(5, 4))
        ax.fill_between(E, flux_class, alpha=1, color='darkblue' , label=r'$\Phi_{\mathrm{classic}}$')
        ax.fill_between(E, flux_bound, alpha=1, color='darkorange', label=r'$\Phi_{\mathrm{bounded}}$')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(1e2, 1e8)
        ax.set_ylim(1e-14, 1e-8)
        #ax.set_title(fr'$\gamma={gamma}, \Phi_0={phi0:.0e}$')
        ax.set_xlabel('Energy [GeV]')
        ax.set_ylabel(r'Flux $\Phi(E)$ [GeV$^{-1}$ cm$^{-2}$ s$^{-1}$]')
        ax.grid(True, ls="--", alpha=0.5)
        ax.legend()

        # Filenames
        filename_base = f"flux_gamma{gamma}_phi{phi0:.0e}"
        png_path = os.path.join(out_dir, f"{filename_base}.png")
        pdf_path = os.path.join(out_dir, f"{filename_base}.pdf")

        # Save files
        plt.savefig(png_path, dpi=300)
        plt.savefig(pdf_path)
        plt.close(fig)

print(f"Saved plots in '{out_dir}' folder.")
