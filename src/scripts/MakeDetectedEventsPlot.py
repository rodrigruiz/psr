#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py
from h5py import File
from astropy.table import Table
from astropy.coordinates import SkyCoord, AltAz
import astropy.units as u
from km3astro.frame import get_location


# ---------------------------
# Placeholder for missing imports
# ---------------------------
def load_aeff(file_path, flavor):
    with h5py.File(file_path, 'r') as f:

        aeff_1d = f[flavor]['aeff_1d'][()]
        aeff_2d = f[flavor]['aeff_2d'][()]
        e_bins = f['energy_bins'][()]
        costh_bins = f['cos_theta_bins'][()]
    return aeff_1d, aeff_2d, e_bins, costh_bins

def powerLawFlux(E, Phi0, gamma, E0=1e5):
    """
    Returns differential flux at energy E [GeV], normalized at E0 [GeV].
    Phi0 in units of GeV^-1 cm^-2 s^-1 sr^-1.
    """
    return Phi0 * (E / E0) ** (-gamma)

import numpy as np

def bounded_power_law_flux(E, Phi0, gamma, E0=1e5, E_min=None, E_max=None):
    """
    Power-law differential flux with hard cutoffs.

    Parameters
    ----------
    E : array_like
        Energy(s) in GeV.
    Phi0 : float
        Normalization at E0 (units GeV^-1 cm^-2 s^-1 sr^-1).
    gamma : float
        Spectral index (flux ~ (E/E0)^(-gamma)).
    E0 : float
        Reference energy in GeV where Phi0 is defined.
    E_min, E_max : float or None
        If provided, flux is zero for E < E_min or E > E_max.

    Returns
    -------
    flux : ndarray
        Differential flux evaluated at E.
    """
    E = np.asarray(E, dtype=float)
    flux = Phi0 * (E / E0) ** (-gamma)

    if E_min is not None:
        flux = np.where(E < E_min, 0.0, flux)
    if E_max is not None:
        flux = np.where(E > E_max, 0.0, flux)

    return flux


# ---------------------------
# Plotting function
# ---------------------------
def plot_event_energies(df_events, energy_bins, total_energies, detectorname, phi0, gamma, plot_config=None):
    if plot_config is None:
        plot_config = {
            "flavors": {
                "nue": {"color": "tab:blue", "label": "νe", "linewidth": 1.2},
                "anue": {"color": "tab:orange", "label": "ν̄e", "linewidth": 1.2},
                "numu": {"color": "tab:green", "label": "νμ", "linewidth": 1.2},
                "anumu": {"color": "tab:red", "label": "ν̄μ", "linewidth": 1.2},
                "nutau": {"color": "tab:purple", "label": "ντ", "linewidth": 1.2},
                "anutau": {"color": "tab:brown", "label": "ν̄τ", "linewidth": 1.2},
            },
            "shower_fill": {"color": "skyblue", "alpha": 0.4, "label": "Total showers"},
            "track_fill": {"color": "salmon", "alpha": 0.4, "label": "Total tracks"},
            "total_line": {"color": "black", "linewidth": 2.0, "label": "Total flux"}
        }

    energy_vals = 0.5 * (energy_bins[:-1] + energy_bins[1:])

    fig, ax1 = plt.subplots(figsize=(7, 4))

    def bin_counts(energies):
        return np.histogram(energies, bins=energy_bins)[0]

    shower_energies = df_events[df_events["event_type"] == "shower"]["energy"]
    track_energies = df_events[df_events["event_type"] == "track"]["energy"]

    n_showers = len(shower_energies)
    n_tracks = len(track_energies)

    shower_counts = bin_counts(shower_energies)
    track_counts = bin_counts(track_energies)

    ax1.fill_between(
        energy_vals, track_counts, step='mid',
        color=plot_config["track_fill"]["color"],
        alpha=plot_config["track_fill"]["alpha"],
        label=f"{plot_config['track_fill']['label']} ({n_tracks} events)"
    )

    ax1.fill_between(
        energy_vals, shower_counts, step='mid',
        color=plot_config["shower_fill"]["color"],
        alpha=plot_config["shower_fill"]["alpha"],
        label=f"{plot_config['shower_fill']['label']} ({n_showers} events)"
    )

    for flavor in df_events["flavor"].unique():
        energies = df_events[df_events["flavor"] == flavor]["energy"]
        n_events_flavor = len(energies)
        style = plot_config["flavors"].get(flavor, {})
        ax1.hist(
            energies, bins=energy_bins, histtype='step',
            linewidth=style.get("linewidth", 1.0),
            color=style.get("color", "gray"),
            label=f"{style.get('label', flavor)} ({n_events_flavor} events)",
            linestyle=style.get("linestyle", "-")
        )

    n_total = len(total_energies)
    ax1.hist(
        total_energies, bins=energy_bins, histtype='step',
        linewidth=plot_config["total_line"]["linewidth"],
        color=plot_config["total_line"]["color"],
        label=f"{plot_config['total_line']['label']} ({n_total} events)",
        linestyle=plot_config["total_line"].get("linestyle", "-")
    )

    # --- New: Add secondary y-axis for flux ---
    ax2 = ax1.twinx()
    #flux_vals = powerLawFlux(energy_vals, phi0, gamma)  # in physical units
    flux_vals = bounded_power_law_flux(energy_vals,phi0,gamma,E_min=5e4,E_max=5e5)
    ax2.plot(
        energy_vals, flux_vals,
        color='forestgreen', linestyle='--', linewidth=2,
        label=f"Input PowerLaw Flux (γ={gamma})"
    )
    ax2.set_ylabel(r"Flux [GeV$^{-1}$ cm$^{-2}$ s$^{-1}$ sr$^{-1}$]")
    ax2.set_yscale('log')

    # --- Formatting ---
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel("Energy [GeV]")
    ax1.set_ylabel("Number of events")
    ax1.set_title(f"Simulated Detected Events & Input Flux ({detectorname}, phi0={phi0}, gamma={gamma})")
    ax1.grid(True, which='both', ls='--', alpha=0.5)

    # Combine legends from both axes
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='best', fontsize='small')

    plt.tight_layout()
    plt.savefig(f"detected_events/detected_energies_with_fluxaxis_phi0{phi0}_{detectorname}_gamma{gamma}_total{n_total}.pdf")
    plt.close(fig)


# ---------------------------
# Simulation helpers
# ---------------------------
def simulate_events_with_2d_aeff(phi0, energy_bins, aeff_2d, theta_bin_edges,
                                 source_coord, times, location,
                                 gamma=1.5, T=1e6, Omega=2*np.pi, n_events=None,
                                 E0=1e5, plot=False):
    from astropy.time import Time

    E_centers = 0.5 * (energy_bins[:-1] + energy_bins[1:])
    dE = np.diff(energy_bins)

    #flux = powerLawFlux(E_centers, phi0, gamma)
    flux = bounded_power_law_flux(E_centers, phi0, gamma, E_min=5e4, E_max = 5e5)
    altaz_frame = AltAz(obstime=Time(times), location=location)
    source_altaz = source_coord.transform_to(altaz_frame)
    theta_vals = (90 * u.deg - source_altaz.alt).to(u.rad).value

    weights = np.zeros_like(E_centers)
    dt = T / len(theta_vals)

    for theta in theta_vals:
        theta_idx = np.digitize(theta, theta_bin_edges) - 1
        theta_idx = np.clip(theta_idx, 0, aeff_2d.shape[1] - 1)
        aeff_theta = aeff_2d[:, theta_idx]
        weights += flux * aeff_theta * dE * dt * Omega

    MAX_EVENTS = int(1e8)
    if n_events is None:
        expected_total = np.sum(weights)
        n_events = np.random.poisson(expected_total)
        if n_events > MAX_EVENTS:
            print(f"Warning: n_events={n_events} is very large; truncating to {MAX_EVENTS}.")
            n_events = MAX_EVENTS

    if n_events == 0:
        if plot:
            print("No events to simulate (n_events=0).")
        return np.array([]), 0

    p_bins = weights / np.sum(weights)
    sampled_bins = np.random.choice(len(p_bins), size=n_events, p=p_bins)

    log_E_low = np.log10(energy_bins[sampled_bins])
    log_E_high = np.log10(energy_bins[sampled_bins + 1])
    simulated_logE = np.random.uniform(log_E_low, log_E_high)
    simulated_energies = 10**simulated_logE

    if plot:
        plt.figure(figsize=(8, 5))
        plt.hist(simulated_energies, bins=energy_bins, histtype='step', color='blue', lw=1.5)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel("Energy [GeV]")
        plt.ylabel("Number of events")
        plt.title(f"Simulated Detected Event Energies\nPhi0={phi0:.1e}, gamma={gamma}")
        plt.grid(True, which='both', ls='--', alpha=0.5)
        plt.tight_layout()
        plt.show()

    return simulated_energies, n_events


def simulate_events_all_flavors(phi0, energy_bins, theta_bin_edges, aeff_dict,
                                source_coord, times, location, detectorname,
                                gamma=1.5, days=90, Omega=2*np.pi,
                                n_events=None, E0=1e5, plot=False, plot_config=None):
    total_events = 0
    total_energies = []
    event_records = []
    shower_energies = []
    track_energies = []

    T = days * 86400

    FLAVOR_EVENTTYPE_MAP = {
        "nue": "shower", "anue": "shower",
        "numu": "track", "anumu": "track",
        "nutau": "shower", "anutau": "shower"
    }

    n_flavors = len(aeff_dict)
    phi0_per_flavor = phi0 / n_flavors

    for flavor, aeff_2d in aeff_dict.items():
        simulated_energies, n_evt = simulate_events_with_2d_aeff(
            phi0=phi0_per_flavor, energy_bins=energy_bins,
            aeff_2d=aeff_2d, theta_bin_edges=theta_bin_edges,
            source_coord=source_coord, times=times, location=location,
            gamma=gamma, T=T, Omega=Omega, E0=E0, plot=False
        )

        total_events += n_evt
        total_energies.extend(simulated_energies)

        event_type = FLAVOR_EVENTTYPE_MAP.get(flavor, "shower")
        if event_type == "shower":
            shower_energies.extend(simulated_energies)
        elif event_type == "track":
            track_energies.extend(simulated_energies)

        records = [{"energy": energy, "flavor": flavor, "event_type": event_type}
                   for energy in simulated_energies]
        event_records.extend(records)

    df_events = pd.DataFrame(event_records)
    #em_shower_events = df_events[df_events["flavor"].isin(["nue", "anue"])].shape[0]

    #print("em_shower_events:", em_shower_events)
    print("df_events:", df_events)

    if plot and len(total_energies) > 0:
        plot_event_energies(df_events, energy_bins, total_energies,
                            detectorname, phi0, gamma, plot_config)

    return (np.array(total_energies), total_events,
            df_events, np.array(shower_energies), np.array(track_energies))


def run_simulation_for_detectors(detectorname, phi0_list, gamma_list,
                                 source_file, plot_config=None):
    if detectorname.lower() == "orca":
        aeff_file = "aeff_summary_v9_flavors_orca.hdf5"
        input_file = "/home/wecapstor3/capn/capn107h/combined_eventlists/orca_KM3NeT_00000148_0001_combined_1008510events_mc.hdf5"
    elif detectorname.lower() == "arca":
        aeff_file = "aeff_summary_v9_flavors_arca.hdf5"
        input_file = "/home/wecapstor3/capn/capn107h/combined_eventlists/arca_KM3NeT_00000133_0001_combined_30929077events_data.hdf5"
    else:
        raise ValueError(f"Detector name '{detectorname}' not recognized.")

    aeff_dict = {}
    for flavor in ['nue', 'anue', 'numu', 'anumu', 'nutau', 'anutau']:
        _, aeff_2d, e_bins, cos_theta_bins = load_aeff(aeff_file, flavor)
        aeff_dict[flavor] = aeff_2d

    with h5py.File(source_file, 'r') as f:
        ra = f['skycoord/ra'][()]
        dec = f['skycoord/dec'][()]
        source_skycoord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

    with File(input_file, 'r') as h5_file:
        eventList = Table.read(h5_file)
        time = eventList['time']

    detector_location = get_location(detectorname)
    observation_times = np.linspace(time.min(), time.max(), 1000)

    energy_bins = e_bins
    theta_bin_edges = np.arccos(cos_theta_bins[::-1])

    for phi0 in phi0_list:
        for gamma in gamma_list:
            energies, n_total, event_df, _, _ = simulate_events_all_flavors(
                phi0=phi0, energy_bins=energy_bins, theta_bin_edges=theta_bin_edges,
                aeff_dict=aeff_dict, source_coord=source_skycoord, times=observation_times,
                location=detector_location, detectorname=detectorname, days=90,
                gamma=gamma, plot=False, plot_config=plot_config
            )

            print(f"[{detectorname} - phi0={phi0}, gamma={gamma}] Total events: {n_total}")
            #
            #print(f"[{detectorname} - phi0={phi0}, gamma={gamma}] EM shower events: {n_em_showers}")

            plot_event_energies(event_df, energy_bins, energies,
                                detectorname, phi0, gamma, plot_config)


# ---------------------------
# Main script entry
# ---------------------------
if __name__ == "__main__":
    plot_config = {
        "flavors": {
            "nue":    {"color": "darkblue", "label": r"$\nu_e$", "linestyle": "-", "linewidth": 1.5},
            "anue":   {"color": "darkblue", "label": r"$\bar\nu_e$", "linestyle": "dotted", "linewidth": 2},
            "numu":   {"color": "firebrick", "label": r"$\nu_\mu$", "linestyle": "-", "linewidth": 1.5},
            "anumu":  {"color": "firebrick", "label": r"$\bar\nu_\mu$", "linestyle": "dotted", "linewidth": 2},
            "nutau":  {"color": "gold", "label": r"$\nu_\tau$", "linestyle": "-", "linewidth": 1.5},
            "anutau": {"color": "gold", "label": r"$\bar\nu_\tau$", "linestyle": "dotted", "linewidth": 2},
        },
        "shower_fill": {"color": "skyblue", "alpha": 0.8, "label": "Total showers"},
        "track_fill": {"color": "salmon", "alpha": 0.4, "label": "Total tracks"},
        "total_line": {"color": "black", "linewidth": 2, "label": "Total flux", "linestyle": "-"}
    }

    phi0_values = [1e-14, 1e-13, 1e-12, 1e-11, 1e-10]
    gamma_values = [1.0, 2.0]
    source_file = "/home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_analysis/inputs/sources/Vela_X-1.h5"

    #run_simulation_for_detectors("orca", phi0_values, gamma_values,source_file=source_file, plot_config=plot_config)

    run_simulation_for_detectors("arca", phi0_values, gamma_values,
                                 source_file=source_file,
                                 plot_config=plot_config)
