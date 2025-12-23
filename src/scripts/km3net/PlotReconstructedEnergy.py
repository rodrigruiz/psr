import numpy as np
import h5py
import matplotlib.pyplot as plt

# ----------------------------------------------------------
# 1) Load response from HDF5 file
# ----------------------------------------------------------
with h5py.File("/home/hpc/capn/capn107h/software/psr/src/scripts/arca_energy_response_histogram_arca_track_pdgid14_-14.hdf5", "r") as f:
    resp = f["hist2d"][:]        # shape (n_true_bins, n_reco_bins)
    bins_true = f["bins_true"][:]
    bins_reco = f["bins_reco"][:]

# ----------------------------------------------------------
# 2) Sampling function
# ----------------------------------------------------------
def sample_reco_energy_array(true_energies, resp, bins_true, bins_reco):
    """
    Sample reconstructed energies for an array of true energies
    based on conditional PDF P(E_reco | E_true).
    """
    true_energies = np.asarray(true_energies)
    reco_samples = np.full(true_energies.shape, np.nan, dtype=float)

    for i, E_true in enumerate(true_energies):
        true_bin_idx = np.digitize(E_true, bins_true) - 1

        if 0 <= true_bin_idx < resp.shape[0]:
            pdf = resp[true_bin_idx]
            if np.any(pdf):  # skip empty rows
                reco_bin_idx = np.random.choice(len(pdf), p=pdf)

                # Log-uniform sample inside reco bin
                lo, hi = bins_reco[reco_bin_idx], bins_reco[reco_bin_idx+1]
                u = np.random.rand()
                reco_samples[i] = lo * (hi/lo) ** u

    failed_mask = np.isnan(reco_samples)
    return reco_samples, failed_mask

# ----------------------------------------------------------
# 3) Generate true energies for testing
# ----------------------------------------------------------
rng = np.random.default_rng(0)
true_energies = 10**rng.uniform(2, 7, size=20000)  # log-uniform 1e2–1e7 GeV

# ----------------------------------------------------------
# 4) Sample reco energies
# ----------------------------------------------------------
reco_samples, failed_mask = sample_reco_energy_array(true_energies, resp, bins_true, bins_reco)

true_valid = true_energies[~failed_mask]
reco_valid = reco_samples[~failed_mask]

# ----------------------------------------------------------
# 5) Plots
# ----------------------------------------------------------

# --- Smearing hexbin ---
plt.figure(figsize=(7,6))
plt.hexbin(true_valid, reco_valid, bins='log', gridsize=100, cmap='inferno',
           extent=[bins_true[0], bins_true[-1], bins_reco[0], bins_reco[-1]])
plt.plot([1e2, 1e7], [1e2, 1e7], 'w--', lw=2, label="E_reco = E_true")
plt.xscale('log')
plt.yscale('log')
plt.xlabel("True Energy $E_{true}$ [GeV]")
plt.ylabel("Reco Energy $E_{reco}$ [GeV]")
plt.title("Response Smearing: True → Reco")
plt.colorbar(label="Counts (log)")
plt.legend()
plt.tight_layout()
plt.savefig("ResponseSmearingTest.png")

# --- Distributions ---
plt.figure(figsize=(7,4))
plt.hist(true_valid, bins=np.logspace(2,7,80), histtype='step', lw=2, label="True Energies")
plt.hist(reco_valid, bins=np.logspace(2,7,80), histtype='step', lw=2, label="Reco Samples")
plt.xscale('log')
plt.xlabel("Energy [GeV]")
plt.ylabel("Counts")
plt.title("Distribution of True vs Reco Energies")
plt.legend()
plt.tight_layout()
plt.savefig("ErecoEtrueHistogramTest.png")

# --- Ratio E_reco / E_true ---
plt.figure(figsize=(7,4))
ratio = reco_valid / true_valid
plt.hist(ratio, bins=np.logspace(-2,2,100), histtype='stepfilled', alpha=0.6, color='C0')
plt.xscale('log')
plt.xlabel("E_reco / E_true")
plt.ylabel("Counts")
plt.title("Energy Response Ratio Distribution")
plt.tight_layout()
plt.savefig("ResponseRatioTest.png")
