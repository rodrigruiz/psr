import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


# -------------------------------------------------
# Load from HDF5
# -------------------------------------------------
def load_aeff_from_hdf5(hdf5_file):
    """
    Loads precomputed 2D and 1D effective areas from HDF5 file.
    Returns:
        results: {flavor: {"aeff_2d": array, "aeff_1d": array}}
        binning: {"energy_bins": array, "cos_theta_bins": array}
    """
    results = {}
    with h5py.File(hdf5_file, "r") as f:
        energy_bins = f["energy_bins"][:]
        cos_theta_bins = f["cos_theta_bins"][:]
        dcos = np.diff(cos_theta_bins)

        for flavor in f.keys():
            if flavor in ["energy_bins", "cos_theta_bins"]:
                continue
            aeff_2d = f[flavor]["aeff_2d"][:]
            aeff_1d = np.sum(aeff_2d * dcos[np.newaxis, :] * 0.5, axis=1)
            results[flavor] = {"aeff_2d": aeff_2d, "aeff_1d": aeff_1d}

    binning = {
        "energy_bins": energy_bins,
        "cos_theta_bins": cos_theta_bins
    }
    return results, binning


# -------------------------------------------------
# 2D Plotting
# -------------------------------------------------
def plot_aeff_2D_grid(results, binning, detector_name, flavors=None,
                      flavor_titles=None, vmin=1e-6, vmax=1e4):
    """
    Plots a 2×2 grid of 2D effective area for selected flavors.
    """
    if flavors is None:
        flavors = ["nue", "anue", "numu", "anumu"]
    if flavor_titles is None:
        flavor_titles = {
            "nue": r"$\nu_e$",
            "anue": r"$\bar{\nu}_e$",
            "numu": r"$\nu_\mu$",
            "anumu": r"$\bar{\nu}_\mu$"
        }

    energy_bins = binning["energy_bins"]
    cos_theta_bins = binning["cos_theta_bins"]
    E_centers = 0.5 * (energy_bins[:-1] + energy_bins[1:])
    cos_centers = 0.5 * (cos_theta_bins[:-1] + cos_theta_bins[1:])

    fig, axes = plt.subplots(2, 2, figsize=(7, 5), sharex=True, sharey=True)
    axes = axes.flatten()

    for i, flavor in enumerate(flavors):
        aeff_2d = results[flavor]["aeff_2d"]

        mesh = axes[i].pcolormesh(
            E_centers, cos_centers, aeff_2d.T,
            shading="auto",
            cmap="viridis",
            norm=LogNorm(vmin=vmin, vmax=vmax)
        )

        axes[i].set_title(flavor_titles[flavor])
        axes[i].set_xscale("log")
        axes[i].set_xlim([energy_bins[0], energy_bins[-1]])
        axes[i].set_ylim([-1, 1])

        if i % 2 == 0:
            axes[i].set_ylabel(r"$\cos(\theta)$")
        if i >= 2:
            axes[i].set_xlabel("Energy [GeV]")

    # Colorbar
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    fig.colorbar(mesh, cax=cbar_ax, label=r"A$_{\mathrm{eff}}$ [m$^2$]")

    fig.suptitle(f"2D Effective Area - {detector_name}") #, fontsize=16)
    #plt.tight_layout(rect=[0, 0, 0.9, 0.95])
    plt.grid(alpha=0.5)
    plt.savefig(f"aeff_2D_{detector_name.lower()}.pdf", dpi=200)
    plt.show()


# -------------------------------------------------
# 1D Plotting
# -------------------------------------------------
def plot_aeff_1D(results, binning, plot_config, detector_name):
    """
    Plots 1D integrated effective area vs energy for given plot_config.

    plot_config: dict of {
        flavor: {
            "label": str,
            "color": str,
            "linestyle": str,
            "alpha": float,
            "linewidth": float
        }
    }
    """
    energy_bins = binning["energy_bins"]
    E_centers = 0.5 * (energy_bins[:-1] + energy_bins[1:])

    plt.figure(figsize=(7, 5))
    for flavor, cfg in plot_config.items():
        plt.plot(
            E_centers,
            results[flavor]["aeff_1d"],
            label=cfg.get("label", flavor),
            color=cfg.get("color", "black"),
            linestyle=cfg.get("linestyle", "-"),
            alpha=cfg.get("alpha", 1.0),
            linewidth=cfg.get("linewidth", 1.5),
            drawstyle="steps-mid"
        )

    plt.xlabel("Energy [GeV]")
    plt.ylabel(r"Integrated $A_{\mathrm{eff}}$ [m$^2$]")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid(True, which="both")
    plt.legend()
    plt.title(f"1D Effective Area - {detector_name}")
    plt.tight_layout()
    plt.savefig(f"aeff_1D_{detector_name.lower()}.pdf", dpi=200)
    plt.show()


# -------------------------------------------------
# Example usage
# -------------------------------------------------
if __name__ == "__main__":
    # Files
    files = {
        "ARCA": "/home/hpc/capn/capn107h/software/psr/src/scripts/aeff_summary_v9_flavors_arca.hdf5",
        "ORCA": "/home/hpc/capn/capn107h/software/psr/src/scripts/aeff_summary_v9_flavors_orca.hdf5"
    }


    # Loop through detectors
    for det_name, file_path in files.items():
        results, binning = load_aeff_from_hdf5(file_path)

        # Plot 2D grids
        plot_aeff_2D_grid(results, binning, det_name)

        # Example 1D plot config
        plot_config = {
            "nue":   {"label": r"$\nu_e$",       "color": "blue",  "linestyle": "-",  "alpha": 0.9, "linewidth": 2},
            "anue":  {"label": r"$\bar{\nu}_e$", "color": "blue",  "linestyle": "--", "alpha": 0.9, "linewidth": 2},
            "numu":  {"label": r"$\nu_\mu$",     "color": "red",   "linestyle": "-",  "alpha": 0.9, "linewidth": 2},
            "anumu": {"label": r"$\bar{\nu}_\mu$","color": "red",  "linestyle": "--", "alpha": 0.9, "linewidth": 2}
        }

        # Plot 1D curves
        plot_aeff_1D(results, binning, plot_config, det_name)
