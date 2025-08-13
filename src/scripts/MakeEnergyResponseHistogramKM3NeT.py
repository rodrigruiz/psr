import numpy as np
from h5py import File
from astropy.table import Table
import matplotlib.pyplot as plt
import os

def get_pdg_label(pdgid):
    if set(pdgid) == {14, -14}:
        return "track"
    elif set(pdgid) == {12, -12, 16, -16}:
        return "shower"
    else:
        return "custom"

def create_energy_histogram(input_file, output_file_prefix, pdgid=[14], bins_true=None, bins_reco=None):
    """Create and save a 2D histogram of true vs. reco energies for specified PDG IDs."""
    with File(input_file, 'r') as h5_file:
        EventList = Table.read(h5_file)
        pdgmask = np.isin(EventList['pdg_id'], pdgid)
        E_mc = EventList[pdgmask]['energy_mc']
        E_reco = EventList[pdgmask]['energy']
        n_gen = EventList[pdgmask]['num_gen_events']
        weights = EventList[pdgmask]['normalized_weight']

    valid_mask = (~np.isnan(E_mc) & ~np.isnan(E_reco) & (E_mc > 0) & (E_reco > 0))
    E_mc_clean = E_mc[valid_mask]
    E_reco_clean = E_reco[valid_mask]

    if bins_true is None:
        bins_true = np.logspace(np.log10(E_mc_clean.min()), np.log10(E_mc_clean.max()), 100)
    if bins_reco is None:
        bins_reco = np.logspace(np.log10(E_reco_clean.min()), np.log10(E_reco_clean.max()), 100)

    hist2d, xedges, yedges = np.histogram2d(E_mc_clean, E_reco_clean, bins=[bins_true, bins_reco])#, weights=weights) #, density = True)

    pdg_str = "_".join(str(pid) for pid in pdgid)
    output_file = f"{output_file_prefix}_pdgid{pdg_str}.hdf5"
    with File(output_file, 'w') as f_out:
        f_out.create_dataset("hist2d", data=hist2d)
        f_out.create_dataset("bins_true", data=xedges)
        f_out.create_dataset("bins_reco", data=yedges)

    print(f"Saved histogram to: {output_file}")
    return output_file

def plot_energy_response_histogram(hist_file, detectorname, pdgid, save_dir, sampled_true_energy = None):
    """Load histogram and plot true vs reco energy, save plot to disk."""
    with File(hist_file, 'r') as f:
        hist2d = f['hist2d'][()]
        bins_true = f['bins_true'][()]
        bins_reco = f['bins_reco'][()]

    pdg_label = get_pdg_label(pdgid)
    title = f"{detectorname.upper()} Energy Response ({pdg_label})"
    fname = os.path.join(save_dir, f"{detectorname}_energy_response_{pdg_label}.pdf")

    plt.figure(figsize=(7, 5))
    plt.pcolormesh(bins_true, bins_reco, hist2d.T, norm=plt.matplotlib.colors.LogNorm(), shading='auto')
    plt.plot([bins_true[0], bins_true[-1]], [bins_true[0], bins_true[-1]], 'k--', label=r'$E_{true} = E_{reco}$')

    if sampled_true_energy is not None:
        plt.axvline(sampled_true_energy, color='darkorange', lw=3, alpha=0.5, label=f"Sampled $E_{{true}}$ = {sampled_true_energy:.1e} GeV")

    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(bins_true[0], bins_true[-1])
    plt.ylim(bins_true[0], bins_true[-1])
    plt.xlabel(r"True Energy $E_{true}$ (GeV)")
    plt.ylabel(r"Reconstructed Energy $E_{reco}$ (GeV)")
    plt.title(title)
    plt.colorbar(label='Counts')
    plt.legend()
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    print(f"Saved energy response plot to: {fname}")

def sample_reco_energy(true_energy, hist2d, bins_true, bins_reco):
    """Sample reco energy from histogram given true energy."""
    true_bin_idx = np.digitize(true_energy, bins_true) - 1
    if true_bin_idx < 0 or true_bin_idx >= hist2d.shape[0]:
        raise ValueError(f"True energy {true_energy} outside histogram range.")

    slice_hist = hist2d[true_bin_idx, :]
    pdf = slice_hist / np.sum(slice_hist) if np.sum(slice_hist) > 0 else np.zeros_like(slice_hist)
    bin_centers_reco = 0.5 * (bins_reco[:-1] + bins_reco[1:])
    return pdf, bin_centers_reco

def plot_reco_pdf_from_histogram(hist_file, detectorname, pdgid, true_energy, save_dir):
    """Load histogram, get PDF for true_energy, plot reco PDF (no sampled line), save plot."""
    with File(hist_file, 'r') as f:
        hist2d = f['hist2d'][()]
        bins_true = f['bins_true'][()]
        bins_reco = f['bins_reco'][()]

    pdf, reco_centers = sample_reco_energy(true_energy, hist2d, bins_true, bins_reco)
    if np.sum(pdf) == 0:
        print(f"No reco PDF available for true energy {true_energy:.2e} GeV in file {hist_file}")
        return

    pdg_label = get_pdg_label(pdgid)
    title = f"{detectorname.upper()} Reco Energy PDF ({pdg_label}) at $E_{{true}}$={true_energy:.2e} GeV"
    fname = os.path.join(save_dir, f"{detectorname}_reco_pdf_{pdg_label}_Etrue_{int(true_energy)}.pdf")

    plt.figure(figsize=(8, 4))
    plt.plot(reco_centers, pdf, drawstyle='steps-mid', color="darkorange", label='Reco PDF')
    plt.xscale('log')
    plt.xlabel(r"Reconstructed Energy $E_{reco}$ (GeV)")
    plt.ylabel("Probability Density")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(fname)
    plt.close()
    print(f"Saved reco PDF plot to: {fname}")

def run_all(detectorname, input_file, output_prefix, pdgid_lists, sampled_true_energy, save_dir):
    """Main runner function to create histograms, plot them, and plot reco PDFs."""

    # Create output dir if missing
    os.makedirs(save_dir, exist_ok=True)

    # Set energy bins based on detector
    if detectorname.lower() == "orca":
        bins_true = np.logspace(0, 4, 100)  # 1e0 to 1e4
        bins_reco = bins_true.copy()
    elif detectorname.lower() == "arca":
        bins_true = np.logspace(2, 8, 100)  # 1e2 to 1e8
        bins_reco = bins_true.copy()
    else:
        raise ValueError(f"Unknown detector name: {detectorname}")

    for pdgid in pdgid_lists:
        pdg_label = get_pdg_label(pdgid)
        output_file = create_energy_histogram(input_file, output_prefix + f"_{detectorname}_{pdg_label}", 
                                              pdgid=pdgid, bins_true=bins_true, bins_reco=bins_reco)
        plot_energy_response_histogram(output_file, detectorname, pdgid, save_dir, sampled_true_energy)
        plot_reco_pdf_from_histogram(output_file, detectorname, pdgid, sampled_true_energy, save_dir)

    print("All done.")


def main():
    input_file_arca  = "/home/wecapstor3/capn/capn107h/combined_eventlists/arca_KM3NeT_00000133_0001_combined_4983084events_mc.hdf5"
    input_file_orca  = "/home/wecapstor3/capn/capn107h/combined_eventlists/orca_KM3NeT_00000148_0001_combined_1008510events_mc.hdf5"

    save_dir = "./energy_plots"
    sampled_true_energy_arca = 1e5  # GeV, can be changed to any value within the hist range
    sampled_true_energy_orca = 1e2
    pdgid_shower = [12, -12, 16, -16]
    pdgid_track = [14, -14]

    # For ARCA
    run_all(
        detectorname="arca",
        input_file=input_file_arca,
        output_prefix="arca_energy_response_histogram",
        pdgid_lists=[pdgid_shower, pdgid_track],
        sampled_true_energy=sampled_true_energy_arca,
        save_dir=save_dir
    )

    # For ORCA
    run_all(
        detectorname="orca",
        input_file=input_file_orca,
        output_prefix="orca_energy_response_histogram",
        pdgid_lists=[pdgid_shower, pdgid_track],
        sampled_true_energy=sampled_true_energy_orca,
        save_dir=save_dir
    )


if __name__ == "__main__":
    main()