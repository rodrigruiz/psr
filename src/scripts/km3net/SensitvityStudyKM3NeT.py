#!/usr/bin/env python3
import argparse

import os
import glob
from docopt import docopt

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import h5py 

import astropy.units as u
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.time import Time
from astropy.table import Table, vstack
from astropy.timeseries import TimeSeries

#from plens.PulseModel import MVMD, sinusoid
from stingray import EventList, Lightcurve

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from stingray.pulse.search import epoch_folding_search, z_n_search, search_best_peaks
from epochfolding.stingray_epochfolding import savehdf5, get_testfrequencies
from epochfolding.gtis import loadGTIs, saveGTIs

from km3astro.frame import get_location
import km3io.definitions as kd

import numpy as np

import matplotlib as mpl
#mpl.rcParams['font.size'] = 11
mpl.rcParams['axes.labelsize'] = 12
#mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['legend.fontsize'] = 11
mpl.rcParams['figure.titlesize'] = 14

def MVMD(t, f, phi, kappa, a, baseline=None):
    """Modivied von Mises distribution (MVMD).
        Reference: "Fourier Techniques for Very Long Astrophysical Time-Series Analysis" 
                    Scott M. Ransom et al 2002 AJ 124 1788
                    
    Parameters
    ----------
        t : array_like
            Values to evaluate the MVMD at.
        kappa : float
            Shape parameter giving the width of the function.
        a : float
            Amplitude of the Pulse. Equates to the area of one pulse.
        f : float
            Frequency of the pulse train.
        phi : float
            Phase offset of the pulse train.
    
    Returns
    -------
        y : 1darray
            Evaluated values of the MVMD.
    
    For kappa -> 0 the MVMD converges towards a sinusod.
    For kappa -> infinity the MVMD converges towards a Gaussian. 1/kappa corresponds to sigma^2. Keep in mind, that  
        
    """
    
    y = a * ( np.exp(kappa*np.cos(2*np.pi*f*t+phi))-np.exp(-kappa))/(np.i0(kappa) - np.exp(-kappa))
    
    if baseline is not None:
        y += baseline
        
    return y

def sinusoid(times, frequency, baseline, amplitude, phase):
    return baseline + amplitude * np.sin(2 * np.pi * (frequency * times + phase))

def load_histogram(file_path):
    with h5py.File(file_path, 'r') as f:
        hist2d = f['hist2d'][()]
        bins_true = f['bins_true'][()]
        bins_reco = f['bins_reco'][()]
    return hist2d, bins_true, bins_reco

def powerLawFlux(E, Phi0, gamma, E0=1e5):
    """
    Returns differential flux at energy E [GeV], normalized at E0 [GeV].
    Phi0 in units of GeV^-1 cm^-2 s^-1 sr^-1.
    """
    return Phi0 * (E / E0) ** (-gamma)

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


def softbounded_power_law_flux(E, phi0, gamma, E0, Emin, Emax, delta = 50, factor = 2, xoffset = 0):
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

def load_aeff(file_path, flavor):
    with h5py.File(file_path, 'r') as f:

        aeff_1d = f[flavor]['aeff_1d'][()]
        aeff_2d = f[flavor]['aeff_2d'][()]
        e_bins = f['energy_bins'][()]
        costh_bins = f['cos_theta_bins'][()]
    return aeff_1d, aeff_2d, e_bins, costh_bins


def simulate_events_with_2d_aeff(
    phi0,
    energy_bins,
    aeff_2d,
    theta_bin_edges,
    source_coord,
    times,
    location,
    gamma=1.5,
    T=1e6,  # Total analysis time [s]
    E0=1e5,  # Reference energy for phi0 [GeV]
    Emin=4e7,
    Emax=6e5,
    diffuse=False,  # If True: treat flux as diffuse (per sr); if False: point source
    Omega=None,  # Solid angle in sr (required if diffuse=True)
    plot=False,
    max_events = 1e5,
    fluxtype='classic',
):
    """
    Simulate detected event energies from a power-law source with 2D effective area.

    Parameters:
        phi0 (float): Flux normalization at E0
                      Units depend on `diffuse`:
                        - Point source: GeV⁻¹ cm⁻² s⁻¹
                        - Diffuse: GeV⁻¹ cm⁻² s⁻¹ sr⁻¹
        energy_bins (array): Energy bin edges [GeV]
        aeff_2d (2D array): Effective area [cm²], shape (n_energy_bins-1, n_theta_bins)
        theta_bin_edges (array): Zenith angle bin edges [rad]
        source_coord (SkyCoord): Source sky coordinates
        times (array): Observation times (astropy-compatible)
        location (EarthLocation): Detector location
        gamma (float): Power-law spectral index
        T (float): Total observation time [s] (should match time span of `times`)
        E0 (float): Reference energy for power-law normalization
        diffuse (bool): Whether the flux is diffuse or from a point source
        Omega (float or None): Solid angle [sr] over which diffuse flux applies
                               Required if diffuse=True
        plot (bool): Plot resulting histogram if True

    Returns:
        simulated_energies (np.ndarray): Detected event energies [GeV]
        n_events (int): Number of simulated events
    """
    from astropy.coordinates import AltAz
    from astropy.time import Time
    import numpy as np
    import matplotlib.pyplot as plt

    if diffuse and Omega is None:
        raise ValueError("If diffuse=True, you must provide Omega (solid angle in sr).")

    # Compute energy bin centers and widths
    E_centers = 0.5 * (energy_bins[:-1] + energy_bins[1:])
    dE = np.diff(energy_bins)

    # Power-law flux at bin centers
    if fluxtype  == 'classic':
        flux = powerLawFlux(E_centers, phi0, gamma, E0=E0)
    elif fluxtype == 'bounded':
        flux = bounded_power_law_flux(E_centers, phi0, gamma, E0=E0, E_min=Emin, E_max=Emax)

    #flux = softbounded_power_law_flux(E_centers, phi0, gamma, E0=E0, E_min=E0, E_max=5e5,delta=10000)
    #flux = bounded_power_law_flux(E_centers, phi0, gamma, E0=E0, E_min=Emin, E_max=Emax)

    # Convert source position to AltAz to get zenith angle at each time
    altaz_frame = AltAz(obstime=Time(times), location=location)
    source_altaz = source_coord.transform_to(altaz_frame)
    theta_vals = (90 * u.deg - source_altaz.alt).to(u.rad).value

    # Initialize expected event weights per energy bin
    weights = np.zeros_like(E_centers)

    dt = T / len(theta_vals)  # Assumes uniform time sampling

    for theta in theta_vals:
        theta_idx = np.digitize(theta, theta_bin_edges) - 1
        theta_idx = np.clip(theta_idx, 0, aeff_2d.shape[1] - 1)
        aeff_theta = aeff_2d[:, theta_idx]

        # Main physics formula:
        # For point source: flux × aeff × dE × dt
        # For diffuse source: flux × aeff × dE × dt × Ω
        if diffuse:
            weights += flux * aeff_theta * dE * dt * Omega
        else:
            weights += flux * aeff_theta * dE * dt

    # Total number of expected events
    expected_total = np.sum(weights)
    n_events = np.random.poisson(expected_total)

    if n_events > max_events:
        raise RuntimeError(
            f"Too many events expected ({n_events:.2e} > {max_events:.2e}). "
            "Skipping simulation."
        )

    if n_events == 0:
        if plot:
            print("No events to simulate (n_events = 0).")
        return np.array([]), 0

    # Normalize to probability distribution across energy bins
    p_bins = weights / np.sum(weights)
    sampled_bins = np.random.choice(len(p_bins), size=n_events, p=p_bins)

    # Sample energies within bins (log-uniformly)
    log_E_low = np.log10(energy_bins[sampled_bins])
    log_E_high = np.log10(energy_bins[sampled_bins + 1])
    simulated_logE = np.random.uniform(log_E_low, log_E_high)
    simulated_energies = 10**simulated_logE

    if plot:
        plt.figure(figsize=(8, 5))
        plt.hist(simulated_energies, bins=energy_bins, histtype='step', lw=1.5)
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel("Energy [GeV]")
        plt.ylabel("Number of events")
        title_type = "Diffuse" if diffuse else "Point Source"
        plt.title(f"Simulated Detected Event Energies ({title_type})\nΦ₀={phi0:.1e}, γ={gamma}")
        plt.grid(True, which='both', ls='--', alpha=0.5)
        plt.tight_layout()
        plt.show()

    return simulated_energies, n_events



def plot_event_counts_vs_energy_and_phi0_with_contours(
    phi0_values,
    energy_bins,
    mean_A_eff,
    gamma=1.5,
    T=1e6,
    Omega=2*np.pi,
    contour_levels=None,
    log_counts=True,
    E0=1e5,
    Emin=7e4,
    Emax=6e5,
    plot=False
):
    """
    Plots a 2D heatmap: event counts per energy bin vs. phi0, with optional contour overlays.

    Parameters:
    - phi0_values: array of flux normalization values [cm⁻² s⁻¹ sr⁻¹ GeV⁻¹]
    - energy_bins: array of energy bin edges [GeV]
    - mean_A_eff: array of mean effective areas per energy bin [cm²]
    - gamma: power-law index
    - T: exposure time [s]
    - Omega: solid angle [sr]
    - contour_levels: list of contour levels (in raw counts, not log10)
    - log_counts: whether to plot log10(counts)
    - E0: normalization energy [GeV]
    - plot: bool, if True will generate the plot
    """
    E_centers = np.sqrt(energy_bins[:-1] * energy_bins[1:])
    dE = np.diff(energy_bins)

    counts_matrix = []
    for phi0 in phi0_values:
        flux = powerLawFlux(E_centers, phi0, gamma) # phi0 * (E_centers / E0) ** (-gamma)
        #flux = bounded_power_law_flux(E_centers, phi0, gamma, E_min = Emin, E_max= Emin)
        counts = flux * mean_A_eff * T * dE * Omega
        counts_matrix.append(counts)

    counts_matrix = np.array(counts_matrix)  # shape: (len(phi0_values), len(energy_bins)-1)

    if not plot:
        return counts_matrix

    # For plotting
    logE = np.log10(E_centers)
    logPhi0 = np.log10(phi0_values)
    X, Y = np.meshgrid(logE, logPhi0)

    Z = counts_matrix
    if log_counts:
        Z = np.where(Z > 0, Z, 1e-10)
        Z_plot = np.log10(Z)
    else:
        Z_plot = Z

    # Plot heatmap
    plt.figure(figsize=(10, 6))
    cmap = 'viridis'
    im = plt.imshow(Z_plot, aspect='auto', origin='lower',
                    extent=[logE.min(), logE.max(), logPhi0.min(), logPhi0.max()],
                    cmap=cmap)

    plt.xlabel("log10(Energy [GeV])")
    plt.ylabel("log10(Phi₀ [cm⁻² s⁻¹ sr⁻¹ GeV⁻¹])")
    cbar_label = "log10(Event Counts)" if log_counts else "Event Counts"
    plt.colorbar(im, label=cbar_label)

    # Add contours if requested
    if contour_levels is not None:
        contour_vals = np.log10(contour_levels) if log_counts else contour_levels
        CS = plt.contour(X, Y, Z_plot, levels=contour_vals, colors='white', linewidths=1.5)
        plt.clabel(CS, fmt=lambda v: f"{10**v:.0f}" if log_counts else f"{v:.0f}", inline=True, fontsize=10)

    plt.title("Expected Event Counts per Energy Bin vs. Phi₀")
    plt.grid(True, which='both', ls='--', alpha=0.4)
    plt.tight_layout()
    plt.show()


def plot_event_energies(df_events, energy_bins, total_energies,phi0,gamma,detector, plot_config=None):
    """
    Plot event energy distributions per flavor, total flux, and shaded regions
    for total showers and tracks.

    Parameters:
        df_events (pd.DataFrame): Must contain 'energy', 'flavor', 'event_type'
        energy_bins (array-like): Energy bin edges
        total_energies (array-like): All energies combined (for total flux)
        plot_config (dict, optional): Style configuration for plotting
    """
    if plot_config is None:
        plot_config = {
            "flavors": {
                "nue":    {"color": "darkblue", "label": r"$\nu_e$",   "linestyle" : "-",   "linewidth": 1.5},
                "anue":   {"color": "darkblue",     "label": r"$\bar\nu_e$",  "linestyle" : "dotted",    "linewidth": 2},
                "numu":   {"color": "firebrick",      "label":r"$\nu_\mu$",  "linestyle" : "-",     "linewidth": 1.5},
                "anumu":  {"color": "firebrick",        "label": r"$\bar\nu_\mu$",  "linestyle" : "dotted",    "linewidth": 2},
                "nutau":  {"color": "gold",     "label": r"$\nu_\tau$",   "linestyle" : "-",    "linewidth": 1.5},
                "anutau": {"color": "gold", "label": r"$\bar\nu_\tau$",   "linestyle" : "dotted",   "linewidth": 2},
            },
            "shower_fill": {
                "color": "skyblue",
                "alpha": 0.8,
                "label": "Total showers"
            },
            "track_fill": {
                "color": "salmon",
                "alpha": 0.4,
                "label": "Total tracks"
            },
            "total_line": {
                "color": "black",
                "linewidth": 2,
                "label": "Total flux",
                "linestyle":  "-"
            }
        }

    # Midpoints for fill_between
    energy_vals = 0.5 * (energy_bins[:-1] + energy_bins[1:])

    plt.figure(figsize=(10, 6))

    

    # Count bin entries
    def bin_counts(energies):
        return np.histogram(energies, bins=energy_bins)[0]

    # Fill: shower and track totals
    shower_energies = df_events[df_events["event_type"] == "shower"]["energy"]
    track_energies = df_events[df_events["event_type"] == "track"]["energy"]

    shower_counts = bin_counts(shower_energies)
    track_counts = bin_counts(track_energies)

    plt.fill_between(
        energy_vals,
        track_counts,
        step='mid',
        color=plot_config["track_fill"]["color"],
        alpha=plot_config["track_fill"]["alpha"],
        label=plot_config["track_fill"]["label"]
    )

    plt.fill_between(
        energy_vals,
        shower_counts,
        step='mid',
        color=plot_config["shower_fill"]["color"],
        alpha=plot_config["shower_fill"]["alpha"],
        label=plot_config["shower_fill"]["label"]
    )

    # Plot each flavor
    for flavor in df_events["flavor"].unique():
        energies = df_events[df_events["flavor"] == flavor]["energy"]
        style = plot_config["flavors"].get(flavor, {})
        plt.hist(
            energies,
            bins=energy_bins,
            histtype='step',
            linewidth=style.get("linewidth", 1.0),
            color=style.get("color", "gray"),
            label=style.get("label", flavor),
            linestyle = style.get("linestyle", "-") 
        )

    # Plot total flux
    plt.hist(
        total_energies,
        bins=energy_bins,
        histtype='step',
        linewidth=plot_config["total_line"]["linewidth"],
        color=plot_config["total_line"]["color"],
        label=plot_config["total_line"]["label"],
        linestyle = plot_config["total_line"]["linestyle"]
    )

    # Final plot styling
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Energy [GeV]")
    plt.ylabel("Number of events")
    plt.title("Simulated Detected Events by Flavor & Type")
    plt.grid(True, which='both', ls='--', alpha=0.5)
    plt.legend(loc='best', fontsize='small')
    plt.tight_layout()
    plt.savefig(f"detected_energies_from_simulated_flux_phi0{phi0}_gamma{gamma}_{detector}_smfig_bounded.png")


def simulate_events_all_flavors(
    phi0,
    energy_bins,
    theta_bin_edges,
    aeff_dict,
    source_coord,
    times,
    location,
    detectorname,
    gamma=1.5,
    days=90,
    Omega=2*np.pi,
    E0=1e5,
    Emin=7e4,
    Emax=6e5,
    plot=False,
    diffuse=False,
    plot_config=None,
    fluxtype='classic'
):
    total_events = 0
    total_energies = []
    event_records = []
    shower_energies = []
    track_energies = []

    T = days * 86400  # observation time in seconds

    FLAVOR_EVENTTYPE_MAP = {
        "nue": "shower",
        "anue": "shower",
        "numu": "track",
        "anumu": "track",
        "nutau": "shower",
        "anutau": "shower"
    }

    n_flavors = len(aeff_dict)
    phi0_per_flavor = phi0 / n_flavors

    for flavor, aeff_2d in aeff_dict.items():
        simulated_energies, n_evt = simulate_events_with_2d_aeff(
            phi0=phi0_per_flavor,
            energy_bins=energy_bins,
            aeff_2d=aeff_2d,
            theta_bin_edges=theta_bin_edges,
            source_coord=source_coord,
            times=times,
            location=location,
            gamma=gamma,
            T=T,
            diffuse=False,
            Omega=None,
            E0=E0,
            Emin=Emin,
            Emax=Emax,
            plot=False,
            fluxtype=fluxtype,
        )

        total_events += n_evt
        total_energies.extend(simulated_energies)

        event_type = FLAVOR_EVENTTYPE_MAP.get(flavor, "shower")
        if event_type == "shower":
            shower_energies.extend(simulated_energies)
        elif event_type == "track":
            track_energies.extend(simulated_energies)

        records = [
            {
                "energy": energy,
                "flavor": flavor,
                "event_type": event_type
            }
            for energy in simulated_energies
        ]
        event_records.extend(records)

    df_events = pd.DataFrame(event_records)

    # Count EM shower events (nu_e and anti_nu_e)
    #em_shower_events = df_events[df_events["flavor"].isin(["nue", "anue"])].shape[0]
    print("Overall events: ", len(total_energies))
    print("Track events: ", len(track_energies))
    print("Shower events: ", len(shower_energies))
    #print("em_shower_events: ", em_shower_events)
    #print("df_events: ", df_events)

    if plot and len(total_energies) > 0:
        plot_event_energies(df_events, energy_bins, total_energies,phi0, gamma, detectorname, plot_config)

    return (
        np.array(total_energies),
        total_events,
        #em_shower_events,
        df_events,
        np.array(shower_energies),
        np.array(track_energies)
    )


def simulate_detected_event_energies_from_flux(
    phi0,
    energy_bins,
    mean_A_eff,
    gamma=1.5,
    T=1e6,
    Omega=2*np.pi,
    n_events=None,
    E0=1e5,
    Emin=7e4,
    Emax=6e5,
    plot=False
):
    """
    Simulates detected event energies based on folded flux × effective area, with optional automatic estimation of n_events.
    
    Parameters:
    - phi0: flux normalization [cm⁻² s⁻¹ sr⁻¹ GeV⁻¹]
    - energy_bins: array of energy bin edges [GeV]
    - mean_A_eff: array of mean effective areas per energy bin [cm²]
    - gamma: spectral index
    - T: exposure time [s]
    - Omega: solid angle [sr]
    - n_events: number of events to simulate (if None, sample from Poisson of expected)
    - E0: normalization energy [GeV]
    - plot: bool, if True will plot the simulated energy distribution

    Returns:
    - simulated_energies: array of sampled event energies [GeV]
    - n_events: number of events simulated
    """
    E_centers = 0.5*(energy_bins[:-1] + energy_bins[1:])
    dE = np.diff(energy_bins)

    flux = powerLawFlux(E_centers, phi0, gamma) # phi0 * (E_centers / E0) ** (-gamma)
    #flux = bounded_power_law_flux(E_centers, phi0, gamma, E_min = Emin, E_max= Emin)
    weights = flux * mean_A_eff * T * dE * Omega

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
        plt.figure(figsize=(8,5))
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




from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# Functions to fit to the angular resolution distribution
def logistic_function(E, a, b, c, d,e,f):
    return a / (1.0 + np.exp(-b * (np.log10(E) - c))) + d + e * np.maximum(0, np.log10(E) - f)

def polynomial(E, a, b, c, d, e, f, g, h):
    return a * np.log10(E)**4 + b * np.log10(E)**3 + c * np.log10(E)**2 + d * np.log10(E) + e + f*np.log10(E)**5 +g * np.log10(E)**6 +h*np.log10(E)**7





# Function to read HDF5 file and fit the data to a model
def fit_angular_resolution(file_path, output_dir, energy_low, energy_high, fit_type="logistic",plot=False):
    # Open the HDF5 file
    table = Table.read(file_path, format='hdf5')

    # Access the columns by their names
    energy_bin_centers = table['energy_bin_center']
    angular_resolution = table['median_angular_separation']
    sigma1_width = table['sigma1_width']  # Assuming this is the uncertainty

    valid_indices = np.isfinite(angular_resolution) & np.isfinite(sigma1_width) & np.isfinite(energy_bin_centers)

    if not np.all(valid_indices):
        print("Warning: Some values were NaN or Inf and have been removed from the fit.")
        energy_bin_centers = energy_bin_centers[valid_indices]
        angular_resolution = angular_resolution[valid_indices]
        sigma1_width = sigma1_width[valid_indices]



    if fit_type == "logistic":
        # Initial guesses for the logistic function parameters
        popt, _ = curve_fit(logistic_function, energy_bin_centers, angular_resolution, 
                            p0=[2, -1, 4, 0.1,1.0, 8], sigma=sigma1_width, absolute_sigma=True)
        fit_func = lambda E: logistic_function(E, *popt)
        plot_name = os.path.join(output_dir,f"fitted_angular_resolution_tracklike_smfig_bounded.png")
    elif fit_type == "polymial":
        # Initial guesses for the sinusoidal function parameters
        popt, _ = curve_fit(polynomial, energy_bin_centers, angular_resolution, p0=[0.0,0.0,0.05,-0.3,1.0,0.0,0.0,0.0] ,sigma=sigma1_width, absolute_sigma=True)
        fit_func = lambda E: polynomial(E, *popt)
        plot_name = os.path.join(output_dir,f"fitted_angular_resolution_showerlike_smfig_bounded.png")
    else:
        raise ValueError("Unknown fit type. Choose 'logistic' or 'polymial'.")
    
    if plot == True:
        # Plot the data and the fit
        plt.figure(figsize=(10, 6))
        plt.plot(energy_bin_centers, angular_resolution, color='black', marker='o', linestyle='none', label='Data')
        energy_fit = np.logspace(energy_low, energy_high, 1000)  # Energy values for plotting the fit
        plt.plot(energy_fit, fit_func(energy_fit), color='firebrick', linestyle='-', linewidth=2, label='Fit')
        plt.xscale('log')
        # Add a shaded 1σ band for uncertainties (optional, depends on your data)
        sigma1_band_lower = angular_resolution - sigma1_width  # Example: replace sigma1_values with your actual 1σ uncertainties
        sigma1_band_upper = angular_resolution + sigma1_width
        plt.fill_between(energy_bin_centers, sigma1_band_lower, sigma1_band_upper, color='blue', alpha=0.3, label='1σ Band')

        # Add a shaded 2σ band for uncertainties (optional)
        sigma2_band_lower = angular_resolution - 2*sigma1_width   # Example: replace sigma2_values with your actual 2σ uncertainties
        sigma2_band_upper = angular_resolution + 2*sigma1_width
        plt.fill_between(energy_bin_centers, sigma2_band_lower, sigma2_band_upper, color='blue', alpha=0.1, label='2σ Band')


        plt.xscale('log')
        # Customize labels and title
        plt.xlabel('Energy [GeV]', fontsize=12)
        plt.ylabel('Angular Resolution [°]', fontsize=12)
        plt.title(f'Angular Resolution vs Energy ({fit_type} Fit Type)', fontsize=14)
        # Add gridlines for better readability
        plt.grid(True, which="both", linestyle="--", linewidth=0.7)
        plt.legend()
        plt.savefig(plot_name)

    return fit_func

def sample_reco_energy_array_old(true_energies, hist2d, bins_true, bins_reco):
    """
    Samples reconstructed energies for an array of true energies based on a 2D histogram.

    Parameters:
    - true_energies: array-like, true neutrino energies
    - hist2d: 2D array of shape (n_true_bins, n_reco_bins), probability distribution
    - bins_true: edges of true energy bins
    - bins_reco: edges of reco energy bins

    Returns:
    - reco_samples: sampled reconstructed energies
    - failed_mask: boolean mask for entries that failed to sample (outside bin range or empty PDF)
    """
    true_energies = np.asarray(true_energies)
    reco_samples = np.full_like(true_energies, np.nan, dtype=float)

    for i, E_true in enumerate(true_energies):
        true_bin_idx = np.digitize(E_true, bins_true) - 1

        # Check if within valid range
        if 0 <= true_bin_idx < hist2d.shape[0]:
            slice_hist = hist2d[true_bin_idx, :]
            if np.sum(slice_hist) > 0:
                pdf = slice_hist / np.sum(slice_hist)
                reco_bin_idx = np.random.choice(len(pdf), p=pdf)

                # Sample uniformly in the selected reco bin
                reco_lo = bins_reco[reco_bin_idx]
                reco_hi = bins_reco[reco_bin_idx + 1]
                reco_samples[i] = np.random.uniform(reco_lo, reco_hi)

    failed_mask = np.isnan(reco_samples)
    return reco_samples, failed_mask

def sample_reco_energy_array(true_energies, resp, bins_true, bins_reco):
    """
    Vectorized version of sample_reco_energy:
    Samples reconstructed energies for an array of true energies.

    Parameters
    ----------
    true_energies : array-like
        True neutrino energies.
    resp : 2D array of shape (n_true_bins, n_reco_bins)
        Response matrix: P(E_reco | E_true) histogram (not necessarily normalized).
    bins_true : array
        Bin edges for true energy.
    bins_reco : array
        Bin edges for reco energy.

    Returns
    -------
    reco_samples : array
        Sampled reconstructed energies (NaN if failed).
    failed_mask : boolean array
        Mask of entries that failed to sample (outside bin range or empty PDF).
    """
    import numpy as np

    true_energies = np.asarray(true_energies)
    reco_samples = np.full_like(true_energies, np.nan, dtype=float)

    for i, E_true in enumerate(true_energies):
        # locate true-energy bin
        t_idx = np.digitize(E_true, bins_true) - 1

        if 0 <= t_idx < resp.shape[0]:
            p = resp[t_idx]
            if np.any(p):
                # normalize row to form PDF
                p = p / np.sum(p)

                # choose reco bin
                r_idx = np.random.choice(len(p), p=p)

                # sample uniformly inside that bin (linear scale, like the first function)
                lo, hi = bins_reco[r_idx], bins_reco[r_idx+1]
                reco_samples[i] = np.random.uniform(lo, hi)

    failed_mask = np.isnan(reco_samples)
    return reco_samples, failed_mask


def sample_reco_directions_from_source_skycoord(
    source_coord,
    reco_energies,
    resolution_func,
    detector_location,
    observation_time
):
    """
    Simulate reconstructed directions (zenith, azimuth) around a source SkyCoord.

    Parameters:
    - source_coord: SkyCoord
        The true source celestial coordinate.
    - reco_energies: np.ndarray
        Reconstructed energies for the injected events.
    - resolution_func: callable
        Function that maps energy to angular resolution in degrees.
    - detector_location: EarthLocation
        Location of the detector.
    - observation_time: Time or float
        Observation time (astropy Time or UNIX timestamp).

    Returns:
    - dict with 'zenith', 'azimuth', 'zenith_true', 'azimuth_true', 'separation'
    """

    # Ensure time is Time object
    if isinstance(observation_time, (float, int)):
        observation_time = Time(observation_time, format='unix')

    n_events = len(reco_energies)
    ang_res_deg = resolution_func(reco_energies)  # e.g., 1-40 deg

    # Create duplicated SkyCoord array for source
    source_coords = SkyCoord(
        ra=np.full(n_events, source_coord.ra.degree) * u.deg,
        dec=np.full(n_events, source_coord.dec.degree) * u.deg,
        frame='icrs'
    )

    # Convert angular resolution (in deg) to radians and sample points
    angular_offsets = np.random.normal(loc=0, scale=ang_res_deg, size=n_events) * u.deg
    position_angles = np.random.uniform(0, 360, size=n_events) * u.deg

    # Offset the source in RA/Dec by these angles
    smeared_coords = source_coords.directional_offset_by(position_angles, angular_offsets)

    # Transform to AltAz (local detector frame)
    altaz_frame = AltAz(obstime=observation_time, location=detector_location)
    smeared_altaz = smeared_coords.transform_to(altaz_frame)
    true_altaz = source_coord.transform_to(altaz_frame)

    azimuths = smeared_altaz.az.to(u.deg).value
    zeniths = (90 * u.deg - smeared_altaz.alt).to(u.deg).value

    # Also return true direction
    az_true = true_altaz.az.to(u.deg).value
    zen_true = (90 * u.deg - true_altaz.alt).to(u.deg).value

    return {
        'zenith': zeniths,
        'azimuth': azimuths,
        'zenith_true': np.full(n_events, zen_true),
        'azimuth_true': np.full(n_events, az_true),
        'skycoord': smeared_coords,
        'separation': smeared_coords.separation(source_coord).to(u.deg).value
    }

def get_source_path(source_coord, detector_location, times):
    """
    Compute the true source track (azimuth, zenith) for a series of times.

    Parameters
    ----------
    source_coord : SkyCoord
        True source coordinates (ICRS).
    detector_location : EarthLocation
        Location of the detector.
    times : astropy Time array
        Array of observation times.

    Returns
    -------
    azimuths : np.ndarray
    zeniths : np.ndarray
    """
    altaz_frame = AltAz(obstime=times, location=detector_location)
    true_altaz = source_coord.transform_to(altaz_frame)

    az = true_altaz.az.to(u.deg).value
    zen = (90 * u.deg - true_altaz.alt).to(u.deg).value
    return az, zen


def plot_reco_directions(result, source_coord, detector_location, times, title="Reconstructed Directions from Vela"):
    zen = result['zenith']
    az = result['azimuth']

    # Get source track
    az_path, zen_path = get_source_path(source_coord, detector_location, times)

    fig, ax = plt.subplots(figsize=(6, 6))

    # Reconstructed events
    ax.scatter(az, zen, s=1, alpha=0.8, label='Reconstructed Events', color='darkblue')

    # Source path (dashed line)
    ax.scatter(az_path, zen_path, linestyle="--", color="darkorange", s=0.3, alpha=0.9, label="Sourcepath")
    #ax.scatter(az_path, zen_path, linestyle="--", color="darkorange", s=2, label="Source")
    # Highlight current true position
    # ax.scatter(result['azimuth_true'][0], result['zenith_true'][0], s=30, color="darkorange", label="Sourcepath")

    ax.set_xlim(0, 360)
    ax.set_ylim(0, 180)
    ax.set_xlabel("Azimuth (°)")
    ax.set_ylabel("Zenith (°)")
    ax.set_title(title)
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.savefig("Testplot_RecoDirs_smfig_bounded.png")
    plt.close()



def inject_signal_times(time, energy, zenith, azimuth, reco_energies_to_inject, reco_types_to_inject, bin_time, pulseshape, frequency,
                        baseline, a, phi, method='classic', kappa=None, num_periods=20, plot=False,
                        skycoord=None, detector_location=None, print_summary=True, preview_n=5, track_resolution_func = None, shower_resolution_func = None, max_separation_deg=5.0, reference_start_time = None, source_name = None, detector_name="earth"):
    """
    Injects signal times into an event list and returns the combined times, energies, and directions.
    Optionally prints a summary of event counts and previews a few example events.

    Parameters:
        ...
        print_summary : bool
            Whether to print a summary of the injection process.
        preview_n : int
            Number of randomly selected events to print as a preview.
    """

    if reference_start_time is None:
        reference_start_time = time.min()  # fallback to old behavior

    # Align to the nearest earlier point on the period grid
    period = 1.0 / frequency
    phase_offset = (time.min() - reference_start_time) % period
    local_start_time = time.min() - phase_offset

    n_events_inject = len(reco_energies_to_inject)
    total_events = len(time)
    num_original_events_to_keep = total_events - n_events_inject

    if num_original_events_to_keep < 0:
        raise ValueError("Number of injected events exceeds original events.")

    # --- Injection process ---
    if method == 'classic':
        high_res_bin_time = 1e-1
        high_res_time = np.arange(local_start_time, time.max(), high_res_bin_time)

        if pulseshape == 'mvm':
            counts = MVMD(high_res_time, frequency, phi, kappa, a, baseline=baseline)
        elif pulseshape == 'sine':
            counts = sinusoid(high_res_time, frequency, baseline, a, phi)
        else:
            raise ValueError("Invalid pulseshape")

        lc_original = Lightcurve(high_res_time, counts, dt=high_res_bin_time / 10, skip_checks=True)
        ev = EventList()
        ev.simulate_times(lc_original)
        new_event_times = ev.time
    
    elif method == 'base':
        start_time = time.min()
        base_interval_length = num_periods / frequency
        high_res_time = np.arange(local_start_time, local_start_time + base_interval_length, bin_time)

        if pulseshape == 'mvm':
            counts = MVMD(high_res_time, frequency, phi, kappa, a, baseline=baseline)
        elif pulseshape == 'sine':
            counts = sinusoid(high_res_time, frequency, baseline, a, phi)
        else:
            raise ValueError("Invalid pulseshape")

        lc_base = Lightcurve(high_res_time, counts, dt=bin_time, skip_checks=True)
        ev = EventList()
        ev.simulate_times(lc_base)
        base_event_times = ev.time

        total_time_span = time.max() - time.min()
        num_intervals = int(np.ceil(total_time_span / base_interval_length))
        chosen_base_times = np.random.choice(base_event_times, n_events_inject, replace=True)
        random_shifts = np.random.randint(0, num_intervals, size=n_events_inject) * base_interval_length
        new_event_times = chosen_base_times + random_shifts

    else:
        raise ValueError("Invalid method: choose 'classic' or 'base'")

    # --- Resample injected events if needed ---
    if len(new_event_times) >= n_events_inject:
        indices_to_keep_new = np.random.choice(len(new_event_times), n_events_inject, replace=False)
    else:
        indices_to_keep_new = np.random.choice(len(new_event_times), n_events_inject, replace=True)
    new_event_times = new_event_times[indices_to_keep_new]
    reco_energies_to_inject = np.array(reco_energies_to_inject)[indices_to_keep_new]
    reco_types_to_inject = np.array(reco_types_to_inject)[indices_to_keep_new]
    # --- Select original events to keep ---
    #indices_to_keep_original = np.random.choice(len(time), num_original_events_to_keep, replace=False)
    #remaining_original_times = time[indices_to_keep_original]
    #remaining_original_energies = energy[indices_to_keep_original]
    #remaining_original_zeniths = zenith[indices_to_keep_original]
    #remaining_original_azimuths = azimuth[indices_to_keep_original]

    assert len(reco_types_to_inject) == len(reco_energies_to_inject), "Mismatch between types and energies"

    def resolution_func_per_event(energies):
        # Input: array of energies
        # Output: array of angular resolutions
        resolutions = np.zeros_like(energies)
        mask_track = np.array(reco_types_to_inject) == 'track'
        mask_shower = np.array(reco_types_to_inject) == 'shower'
        
        resolutions[mask_track] = track_resolution_func(np.array(energies)[mask_track])
        resolutions[mask_shower] = shower_resolution_func(np.array(energies)[mask_shower])
        
        return resolutions

    # --- Sample directions for injected events ---
    reco_dirs = sample_reco_directions_from_source_skycoord(
        source_coord=skycoord,
        reco_energies=reco_energies_to_inject,
        #resolution_func=lambda e: 0.5,  # 5 deg resolution for now
        resolution_func = resolution_func_per_event,
        detector_location=detector_location,
        observation_time=Time(new_event_times, format='unix')
    )


    

    plot = True
    if plot == True:
        zen = reco_dirs['zenith']
        az = reco_dirs['azimuth']
        source_time = Time(np.linspace(new_event_times.min(),new_event_times.max(),200), format='unix')
        # Get source track
        az_path, zen_path = get_source_path(skycoord, detector_location, source_time) # Time(new_event_times, format='unix'))

        fig, ax = plt.subplots(figsize=(7, 5))

        # Reconstructed events
        ax.scatter(az_path, zen_path, linestyle="--", color="darkorange", s=3, alpha=1, label="Sourcepath")

        ax.scatter(az, zen, s=3, alpha=0.8, label=f'Reconstructed Events ({len(az)})', color='darkblue')

        # Source path (dashed line)
        #ax.scatter(az_path, zen_path, linestyle="--", color="darkorange", s=2, label="Source")
        # Highlight current true position
        #ax.scatter(reco_dirs['azimuth_true'][0], reco_dirs['zenith_true'][0], s=30, color="darkorange", label="Sourcepath")
        ax.axhline(90,linestyle="--", color="k",label="Horizon")

        ax.set_xlim(0, 360)
        ax.set_ylim(180, 0) # Top:0, Bottom: 180

        ax.set_xlabel("Azimuth (°)")
        ax.set_ylabel("Zenith (°)")
        ax.set_title(f"Reconstructed Directions from {source_name} (seen from {detector_name})")
        ax.legend()
        ax.grid(True)
        plt.tight_layout()
        plt.savefig(f"Testplot_RecoDirs_{source_name}_{detector_name}_smfig_presentation.png")
        plt.close()

    radec = True
    if radec:
        location = detector_location  # EarthLocation object you already have

        # Observation times of your reconstructed events
        event_times = Time(new_event_times, format='unix')

        # Convert zenith -> altitude
        alt = 90*u.deg - reco_dirs['zenith']*u.deg
        az = reco_dirs['azimuth']*u.deg

        # Build AltAz frame for events
        altaz_frame = AltAz(obstime=event_times, location=location)
        events_altaz = SkyCoord(az=az, alt=alt, frame=altaz_frame)

        # Transform to RA/Dec
        events_icrs = events_altaz.icrs
        ra = events_icrs.ra.deg
        dec = events_icrs.dec.deg

        # --- Source path ---
        source_altaz = SkyCoord(az=az_path*u.deg,
                                alt=(90-zen_path)*u.deg,
                                frame=AltAz(obstime=source_time, location=location))
        source_icrs = source_altaz.icrs
        ra_path = source_icrs.ra.deg
        dec_path = source_icrs.dec.deg

        # --- True source position ---
        true_alt = 90*u.deg - reco_dirs['zenith_true'][0]*u.deg
        true_az  = reco_dirs['azimuth_true'][0]*u.deg
        true_altaz = SkyCoord(az=true_az, alt=true_alt,
                            frame=AltAz(obstime=event_times[0], location=location))
        true_icrs = true_altaz.icrs

        fig, ax = plt.subplots(figsize=(7, 5))

        # Reconstructed events
        ax.scatter(ra, dec, s=3, alpha=0.8, color="darkblue",
                label=f"Reconstructed Events ({len(ra)})")

        # Source path
        #ax.plot(ra_path, dec_path, "--", color="darkorange", lw=1, label="Source path")

        # True source position
        src_ra = true_icrs.ra.deg
        src_dec = true_icrs.dec.deg
        ax.scatter(src_ra, src_dec, s=50, color="darkorange", marker="x", label=f"{source_name}")
        
        for radius, ls in zip([5, 15, 30], ['-', '--', ':']):  # degrees
            circ = Circle((src_ra, src_dec), radius, transform=ax.transData,
                        edgecolor="gray", linestyle=ls, linewidth=1,
                        fill=False, alpha=0.7, label=f"{radius}°")
            ax.add_patch(circ)

        # --- Window centered on source ---
        margin = 35  # degrees, to fit 30° ring
        ax.set_xlim(src_ra + margin+10, src_ra - margin-10)  # reversed RA for sky convention (east left)
        ax.set_ylim(src_dec - margin, src_dec + margin)

        ax.set_xlabel("Right Ascension (°)")
        ax.set_ylabel("Declination (°)")
        ax.set_title(f"Reconstructed Directions around {source_name} (RA/Dec)")
        ax.legend()
        ax.grid(True)

        plt.tight_layout()
        plt.savefig(f"Testplot_RecoDirs_RADEC_{source_name}_{detector_name}_zoom_smfig_presentation.png")
        plt.close()
    


    if max_separation_deg is not None:
        sep = np.array(reco_dirs['separation'])  # in degrees already
        #print("Seperations: ", sep)
        keep_mask = sep <= max_separation_deg
        print(f"len(keep_mask): {len(keep_mask)}/len(sep): {len(sep)} for max sep angle {max_separation_deg} degrees")
        # Filter injected events by the mask
        new_event_times = np.array(new_event_times)[keep_mask]
        reco_energies_to_inject = np.array(reco_energies_to_inject)[keep_mask]
        reco_types_to_inject = np.array(reco_types_to_inject)[keep_mask]
        new_event_zeniths = np.array(reco_dirs['zenith'])[keep_mask]
        new_event_azimuths = np.array(reco_dirs['azimuth'])[keep_mask]

        if print_summary:
            print(f"Applied max separation cut: {max_separation_deg}°")
            print(f"Injected events kept after cut: {keep_mask.sum()} / {len(keep_mask)}")

    else:
        new_event_zeniths = reco_dirs['zenith']
        new_event_azimuths = reco_dirs['azimuth']
    
    #new_event_zeniths = reco_dirs['zenith']
    #new_event_azimuths = reco_dirs['azimuth']

    n_injected_kept = len(new_event_times)
    num_original_events_to_keep = total_events - n_injected_kept
    indices_to_keep_original = np.random.choice(len(time), num_original_events_to_keep, replace=False)
    remaining_original_times = time[indices_to_keep_original]
    remaining_original_energies = energy[indices_to_keep_original]
    remaining_original_zeniths = zenith[indices_to_keep_original]
    remaining_original_azimuths = azimuth[indices_to_keep_original]

    # --- Combine and sort ---
    combined_times = np.concatenate((remaining_original_times, new_event_times))
    combined_energies = np.concatenate((remaining_original_energies, reco_energies_to_inject))
    combined_zeniths = np.concatenate((remaining_original_zeniths, new_event_zeniths))
    combined_azimuths = np.concatenate((remaining_original_azimuths, new_event_azimuths))

    sort_indices = np.argsort(combined_times)
    combined_times = combined_times[sort_indices]
    combined_energies = combined_energies[sort_indices]
    combined_zeniths = combined_zeniths[sort_indices]
    combined_azimuths = combined_azimuths[sort_indices]

    # --- Summary Output ---
    if print_summary:
        print("\n--- Signal Injection Summary ---")
        print(f"Total original events:      {total_events}")
        print(f"Simulated (injected) events: {n_injected_kept} ({n_injected_kept / total_events * 100:.2f}%)")
        print(f"Original events kept:        {num_original_events_to_keep} ({num_original_events_to_keep / total_events * 100:.2f}%)")
        print(f"Combined total events:       {len(combined_times)}")

        # Preview a few combined events
        #print(f"\n--- Preview of {preview_n} Random Combined Events ---")
        #preview_indices = np.random.choice(len(combined_times), min(preview_n, len(combined_times)), replace=False)
        #for idx in preview_indices:
        #    print(f"[{idx}] Time: {combined_times[idx]:.2f}, Energy: {combined_energies[idx]:.2f}, "
        #          f"Zenith: {combined_zeniths[idx]:.2f}, Azimuth: {combined_azimuths[idx]:.2f}")

    return combined_times, combined_energies, combined_zeniths, combined_azimuths



def apply_selection_cuts(event_list, detector, energy_min, energy_max, nhits_tr_min,
                         nhits_sh_min, lik_tr_min, lik_sh_min, beta0_tr_max,
                         zenith_min, trackscore_low_max, trackscore_high_min):
    """Apply a set of selection cuts to a KM3NeT event list."""
    original_lenth = len(event_list['time'])
    # Energy cut
    energy_mask = (event_list['energy'].value >= energy_min) & (event_list['energy'].value <= energy_max)
    event_list = event_list[energy_mask]

    # Track/Shower event selection
    track_events = (event_list['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & \
                   (event_list['rec_stage'] == kd.reconstruction.JMUONBEGIN)

    if detector == 'arca':
        shower_events = (event_list['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & \
                        (event_list['rec_stage'] == kd.reconstruction.AASHOWERBEGIN)
        hits_mask_shower = shower_events & (event_list['AASHOWERFIT_NUMBER_OF_HITS'] >= nhits_sh_min)
    elif detector == 'orca':
        shower_events = (event_list['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & \
                        (event_list['rec_stage'] == kd.reconstruction.JSHOWERBEGIN)
        hits_mask_shower = shower_events & (event_list['JGANDALF_NUMBER_OF_HITS'] >= nhits_sh_min)
    else:
        raise ValueError(f"Unsupported detector: {detector}")

    # Hits cut
    hits_mask_track = track_events & (event_list['JGANDALF_NUMBER_OF_HITS'] >= nhits_tr_min)
    event_list = event_list[hits_mask_shower | hits_mask_track]

    # Re-derive events post hits cut
    track_events = (event_list['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & \
                   (event_list['rec_stage'] == kd.reconstruction.JMUONBEGIN)
    shower_events = (event_list['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & \
                    (event_list['rec_stage'] == kd.reconstruction.AASHOWERBEGIN) if detector == 'arca' else \
                    (event_list['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & \
                    (event_list['rec_stage'] == kd.reconstruction.JSHOWERBEGIN)

    # Likelihood cut
    lik_mask_track = track_events & (np.abs(event_list['likelihood']) >= lik_tr_min)
    lik_mask_shower = shower_events & (np.abs(event_list['likelihood']) >= lik_sh_min)
    event_list = event_list[lik_mask_track | lik_mask_shower]

    # Zenith cut
    zenith_min_rad = zenith_min * np.pi / 180.0
    event_list = event_list[event_list['zenith'] > zenith_min_rad]
    #print(event_list['zenith'])

    # Beta0 cut (only applies to tracks)
    track_events = (event_list['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & \
                   (event_list['rec_stage'] == kd.reconstruction.JMUONBEGIN)
    shower_events = (event_list['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & \
                    (event_list['rec_stage'] == kd.reconstruction.AASHOWERBEGIN) if detector == 'arca' else \
                    (event_list['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & \
                    (event_list['rec_stage'] == kd.reconstruction.JSHOWERBEGIN)
    beta0_tr_mask = track_events & (event_list['JGANDALF_BETA0_RAD'] <= beta0_tr_max)
    event_list = event_list[beta0_tr_mask | shower_events]

    # Trackscore cut
    trackscore_mask = (event_list['track_score'] <= trackscore_low_max) | \
                      (event_list['track_score'] >= trackscore_high_min)
    event_list = event_list[trackscore_mask]

    print(f"New Length adter Cuts: {len(event_list['time'])} of {original_lenth} events ({100*len(event_list['time'])/original_lenth:.1f} %)")

    return event_list

def extract_subset_within_gtis(event_list, gti_array, desired_length_days):
    """Extract a time-limited subset from event_list, using GTIs to accumulate real live-time."""
    subset_event_indices = []
    accumulated_live_time = 0.0
    used_gtis = []
    #print("gti_array in subset function: ", gti_array)
    event_times = event_list['time']
    #print("even_times: ", event_times)
    if gti_array is not None:
        for start, stop in gti_array:
            gti_duration = stop - start
            #print("start: ", start)
            #print("stop: ", stop)
            #print("gti_duration: ", gti_duration)
            
            if accumulated_live_time + gti_duration > desired_length_days:
                remaining_time = desired_length_days - accumulated_live_time
                new_stop = start + remaining_time
                in_partial_gti = (event_times >= start) & (event_times < new_stop)
                selected = np.where(in_partial_gti)[0]
                subset_event_indices.extend(selected.tolist())
                used_gtis.append([start, new_stop])
                accumulated_live_time += remaining_time
                break
            else:
                in_gti = (event_times >= start) & (event_times <= stop)
                selected = np.where(in_gti)[0]
                subset_event_indices.extend(selected.tolist())
                used_gtis.append([start, stop])
                accumulated_live_time += gti_duration

            if accumulated_live_time >= desired_length_days:
                break
    else:
        first_time = event_times[0]
        last_time = first_time + desired_length_days
        in_range = (event_times >= first_time) & (event_times < last_time)
        subset_event_indices = np.where(in_range)[0].tolist()
        used_gtis = [[first_time.value, last_time.value]]
        accumulated_live_time = desired_length_days

    subset = event_list[subset_event_indices]
    return subset, accumulated_live_time, np.array(used_gtis)

def load_static_inputs(input_file, source_file, aeff_file,shower_histo_file, track_histo_file,file_path_angres_sh,file_path_angres_tr, detectorname):
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    print("Loading static inputs...")

    # Load effective areas

    aeff_dict = {}
    for flavor in ['nue', 'anue', 'numu', 'anumu', 'nutau', 'anutau']:
        _, aeff_2d, e_bins, cos_theta_bins = load_aeff(aeff_file, flavor)
        aeff_dict[flavor] = aeff_2d

    # Load energy response histograms
    shower_hist2d, shower_bins_true, shower_bins_reco = load_histogram(shower_histo_file)
    track_hist2d, track_bins_true, track_bins_reco = load_histogram(track_histo_file)

    # Load angular resolution fits
    fit_func_tr = fit_angular_resolution(file_path_angres_tr, ".", np.min(np.log10(e_bins)), np.max(np.log10(e_bins)), fit_type="logistic")
    fit_func_sh = fit_angular_resolution(file_path_angres_sh, ".", np.min(np.log10(e_bins)), np.max(np.log10(e_bins)), fit_type="polymial")

    # Load source coordinates
    with h5py.File(source_file, 'r') as f:
        ra = f['skycoord/ra'][()]
        dec = f['skycoord/dec'][()]
        source_skycoord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

    # Load input events
    
    with h5py.File(input_file) as f:
        eventList = Table.read(f)
        times = eventList['time'] #.value.astype(float)
        observation_times = np.linspace(times.min(), times.max(), 1000)

    detector_location = get_location(detectorname)

    print("Static inputs loaded.")

    return {
        'aeff_dict': aeff_dict,
        'e_bins': e_bins,
        'cos_theta_bins': cos_theta_bins,
        'shower_hist2d': shower_hist2d,
        'shower_bins_true': shower_bins_true,
        'shower_bins_reco': shower_bins_reco,
        'track_hist2d': track_hist2d,
        'track_bins_true': track_bins_true,
        'track_bins_reco': track_bins_reco,
        'fit_func_tr': fit_func_tr,
        'fit_func_sh': fit_func_sh,
        'source_skycoord': source_skycoord,
        'eventList': eventList,
        'times': times,
        'detector_location': detector_location,
        'observation_times': observation_times,
        'detectorname' : detectorname
    }

def load_static_inputs(input_file, source_file, aeff_file, shower_histo_file, track_histo_file, file_path_angres_sh, file_path_angres_tr, detector_name):
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    print(f"Loading static inputs for {detector_name}...")

    # Load effective areas
    aeff_dict = {}
    for flavor in ['nue', 'anue', 'numu', 'anumu', 'nutau', 'anutau']:
        _, aeff_2d, e_bins, cos_theta_bins = load_aeff(aeff_file, flavor)
        aeff_dict[flavor] = aeff_2d

    # Load energy response histograms
    shower_hist2d, shower_bins_true, shower_bins_reco = load_histogram(shower_histo_file)
    track_hist2d, track_bins_true, track_bins_reco = load_histogram(track_histo_file)

    # Load angular resolution fits
    if detector_name == 'arca':
        fit_func_tr = fit_angular_resolution(file_path_angres_tr, ".", 2, 8, fit_type="logistic")
        fit_func_sh = fit_angular_resolution(file_path_angres_sh, ".", 2, 8, fit_type="polymial")
    else: 
        fit_func_tr = fit_angular_resolution(file_path_angres_tr, ".", 0, 2, fit_type="polymial")
        fit_func_sh = fit_angular_resolution(file_path_angres_sh, ".", 0, 2, fit_type="polymial")

    # Load source coordinates
    with h5py.File(source_file, 'r') as f:
        ra = f['skycoord/ra'][()]
        dec = f['skycoord/dec'][()]
        source_skycoord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

    # Load input events
    with h5py.File(input_file) as f:
        eventList = Table.read(f)
        times = eventList['time']
        detector_location = get_location(detector_name)
        observation_times = np.linspace(times.min(), times.max(), 1000)

    print(f"Static inputs (fun2) loaded for {detector_name}.")

    return {
        'aeff_dict': aeff_dict,
        'e_bins': e_bins,
        'cos_theta_bins': cos_theta_bins,
        'shower_hist2d': shower_hist2d,
        'shower_bins_true': shower_bins_true,
        'shower_bins_reco': shower_bins_reco,
        'track_hist2d': track_hist2d,
        'track_bins_true': track_bins_true,
        'track_bins_reco': track_bins_reco,
        'fit_func_tr': fit_func_tr,
        'fit_func_sh': fit_func_sh,
        'source_skycoord': source_skycoord,
        'eventList': eventList,
        'times': times,
        'detector_location': detector_location,
        'observation_times': observation_times,
        'detector_name' : detector_name
    }

def run_detector_analysis(phi0, gamma, E0, static, frequency, max_angle, length_days, reference_start_time,Emin=7e4,Emax=6e5,source_name = "Vela", zenith_cut=90, fluxtype="classic"):
    # Same unpacking and simulation logic from run_analysis up to cut stage
    aeff_dict = static['aeff_dict']
    e_bins = static['e_bins']
    cos_theta_bins = static['cos_theta_bins']
    shower_hist2d = static['shower_hist2d']
    shower_bins_true = static['shower_bins_true']
    shower_bins_reco = static['shower_bins_reco']
    track_hist2d = static['track_hist2d']
    track_bins_true = static['track_bins_true']
    track_bins_reco = static['track_bins_reco']
    fit_func_tr = static['fit_func_tr']
    fit_func_sh = static['fit_func_sh']
    source_skycoord = static['source_skycoord']
    eventList = static['eventList']
    observation_times = static['observation_times']
    detector_location = static['detector_location']
    detector_name = static['detector_name']

    energies, _, _, shower_energies, track_energies = simulate_events_all_flavors(
        phi0=phi0,
        energy_bins=e_bins,
        theta_bin_edges=np.arccos(cos_theta_bins[::-1]),
        aeff_dict=aeff_dict,
        source_coord=source_skycoord,
        times=observation_times,
        location=detector_location,
        detectorname=detector_name,
        days=length_days,
        gamma=gamma,
        plot=False,
        E0=E0,
        Emin=Emin,
        Emax=Emax,
        fluxtype=fluxtype
    )
    print(f"-----> {detector_name} Events Simulated...")

    



    # Check if we actually got any simulated events
    if len(shower_energies) + len(track_energies) == 0:
        print(f"-----> No simulated events for phi0={phi0}, gamma={gamma}. Using original data only.")
        InjectedEventList = eventList.copy()
    else:

        accuracy = 0.5  # in [0.5, 1]
        rng = np.random.default_rng(42)  # for reproducibility
        p_flip = 1.0 - accuracy

        shower_flip_mask = rng.random(len(shower_energies)) < p_flip
        track_flip_mask  = rng.random(len(track_energies))  < p_flip

        def safe_sample_reco_energy_array(energies, hist2d, bins_true, bins_reco):
            if len(energies) == 0:
                return np.array([]), None
            return sample_reco_energy_array(energies, hist2d, bins_true, bins_reco)

        # then replace in your code:
        reco_shower_correct, _ = safe_sample_reco_energy_array(
            shower_energies[~shower_flip_mask], shower_hist2d, shower_bins_true, shower_bins_reco
        )
        reco_shower_mis, _ = safe_sample_reco_energy_array(
            shower_energies[shower_flip_mask],  track_hist2d,  track_bins_true,  track_bins_reco
        )
        reco_track_correct, _ = safe_sample_reco_energy_array(
            track_energies[~track_flip_mask], track_hist2d, track_bins_true, track_bins_reco
        )
        reco_track_mis, _ = safe_sample_reco_energy_array(
            track_energies[track_flip_mask],  shower_hist2d, shower_bins_true, shower_bins_reco
        )


        df_combined = pd.concat([
            # true showers
            pd.DataFrame({
                'reco_energy': reco_shower_correct,
                'true_type':   'shower',
                'type':   'shower',
                'correct':     True
            }),
            pd.DataFrame({
                'reco_energy': reco_shower_mis,
                'true_type':   'shower',
                'type':   'track',
                'correct':     False
            }),
            # true tracks
            pd.DataFrame({
                'reco_energy': reco_track_correct,
                'true_type':   'track',
                'type':   'track',
                'correct':     True
            }),
            pd.DataFrame({
                'reco_energy': reco_track_mis,
                'true_type':   'track',
                'type':   'shower',
                'correct':     False
            }),
        ], ignore_index=True)

        # --- Print summary ---
        n_total   = len(df_combined)
        n_correct = df_combined['correct'].sum()
        n_wrong   = n_total - n_correct

        acc_realized = n_correct / n_total if n_total > 0 else float('nan')

        print("\n===== Misclassification Summary =====")
        print(f"Requested accuracy:  {accuracy:.3f}")
        print(f"Realized accuracy:   {acc_realized:.3f}")
        print(f"Total events:        {n_total}")
        print(f"  Correctly classified:   {n_correct} ({n_correct/n_total:.1%})")
        print(f"  Misclassified:          {n_wrong} ({n_wrong/n_total:.1%})")

        # Breakdown by true type
        for t in ['shower', 'track']:
            df_sub = df_combined[df_combined['true_type'] == t]
            n_sub  = len(df_sub)
            if n_sub == 0:
                continue
            n_sub_corr = df_sub['correct'].sum()
            n_sub_wrong = n_sub - n_sub_corr
            print(f"    {t.capitalize()}s: {n_sub_corr} correct, {n_sub_wrong} misclassified "
                f"({n_sub_corr/n_sub:.1%} correct)")
        print("=====================================\n")


        '''
        reco_shower_energies, _ = sample_reco_energy_array(shower_energies, shower_hist2d, shower_bins_true, shower_bins_reco)
        print("-----> Reconstructed Shower Events Sampled...")
        reco_track_energies, _ = sample_reco_energy_array(track_energies, track_hist2d, track_bins_true, track_bins_reco)
        print("-----> Reconstructed Track Events Sampled...")
        df_combined = pd.concat([
            pd.DataFrame({'reco_energy': reco_shower_energies, 'type': 'shower'}),
            pd.DataFrame({'reco_energy': reco_track_energies, 'type': 'track'})
        ], ignore_index=True)
        '''



        plotdist = False
        if plotdist == True:
            # --- Distributions ---
            true_valid = shower_energies[shower_energies>0]
            reco_valid = reco_shower_energies[reco_shower_energies>0]
            plt.figure(figsize=(7,5))
            if detector_name == 'orca': bins=np.logspace(0,4,50)
            else: bins=np.logspace(2,8,100)
            plt.hist(true_valid, bins=bins, color="darkblue", histtype='step', lw=2, label="True Energies")
            plt.hist(reco_valid, bins=bins, color="darkorange", histtype='step', lw=2, label="Reco Samples")
            plt.xscale('log')
            plt.xlabel("Energy [GeV]")
            plt.ylabel("Counts")
            plt.title(f"Distribution of True vs Reco Energies ({phi0:.2e}, {gamma:.2f}, {detector_name}, shower)")
            plt.grid(which="both",alpha=0.5)
            plt.legend()
            plt.tight_layout()
            plt.savefig(f"EE_HistogramTest_{phi0}_{gamma}_{detector_name}_shower_smfig_bounded.png")

            true_valid = track_energies[track_energies>0]
            reco_valid = reco_track_energies[reco_track_energies>0]
            plt.figure(figsize=(7,5))
            plt.hist(true_valid, bins=bins, color="darkblue", histtype='step', lw=2, label="True Energies")
            plt.hist(reco_valid, bins=bins, color="darkorange", histtype='step', lw=2, label="Reco Samples")
            plt.xscale('log')
            plt.xlabel("Energy [GeV]")
            plt.ylabel("Counts")
            plt.title(f"Distribution of True vs Reco Energies ({phi0:.2e}, {gamma:.2f}, {detector_name}, track)")
            plt.grid(which="both",alpha=0.5)
            plt.legend()
            plt.tight_layout()
            plt.savefig(f"EE_HistogramTest_{phi0}_{gamma}_{detector_name}_track_smfig_bounded.png")

        injected_energies = df_combined['reco_energy']
        injected_types = df_combined['type']

        combined_times, combined_energies, combined_zeniths, combined_azimuths = inject_signal_times(
            time=eventList['time'].value.astype(float),
            energy=eventList['energy'].value.astype(float),
            zenith=eventList['zenith'],
            azimuth=eventList['azimuth'],
            reco_energies_to_inject=injected_energies,
            reco_types_to_inject=injected_types,
            bin_time=1e-5,
            pulseshape='mvm',
            frequency=frequency,
            baseline=0.0,
            a=1.0,
            phi=0.0,
            method='base',
            kappa=5.0,
            skycoord=source_skycoord,
            detector_location=detector_location,
            track_resolution_func=fit_func_tr,
            shower_resolution_func=fit_func_sh,
            max_separation_deg=max_angle,
            reference_start_time=reference_start_time,
            source_name = source_name,
            detector_name=detector_name
        )
        print(f"-----> Signal Injected to {detector_name} data for phi0 {phi0} and gamma {gamma}...")

        InjectedEventList = eventList.copy()
        InjectedEventList['time'] = TimeSeries(time=Time(combined_times, format='unix'))['time']
        InjectedEventList['energy'] = combined_energies
        InjectedEventList['zenith'] = combined_zeniths
        InjectedEventList['azimuth'] = combined_azimuths

    # Apply cuts
    cut_params = {
        'detector': detector_name, 
        'energy_min': 1e4,
        'energy_max': 1e8,
        'nhits_tr_min': 20,
        'nhits_sh_min': 20,
        'lik_tr_min': 50.,
        'lik_sh_min': 50.,
        'beta0_tr_max': 1,
        'zenith_min': zenith_cut,
        'trackscore_low_max': 0.5,
        'trackscore_high_min': 0.5,
    }

    return InjectedEventList
    #return apply_selection_cuts(InjectedEventList, **cut_params)


def main():

    parser = argparse.ArgumentParser(description="Run Sensitivity Study for given source, zenith_cut, opening_angle")
    parser.add_argument("--source", default="vela", help="source name (vela or crab)")
    parser.add_argument("--angle", default=30, help="Opening angle")
    parser.add_argument("--zenith_cut", default=90, help="Maximum zenith angle")
    parser.add_argument("--fluxtype", default="classic", help="Classic or Bounded Power Law Flux (classic or bounded)")
    args = parser.parse_args()

    from datetime import datetime
    datetime_str = datetime.now().strftime('%Y-%m-%d_%H-%M')
    # === User settings ===
    selected_detectors = ['arca'] #, 'orca', 'combined']
    
    #source = "vela"
    #angle = 30
    #zenith_cut = 90 

    source = str(args.source)
    angle = float(args.angle)
    zenith_cut = float(args.zenith_cut) 
    fluxtype = str(args.fluxtype)

    print(f"source : {source}")
    print(f"angle : {angle}")
    print(f"zenith_cut : {zenith_cut}")

    # Optional GTI paths (set to None if not needed)
    #gti_path_arca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_12858events_gtis.pkl"   # or None
    #gti_path_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_58211events_gtis.pkl"   # or None

    #gti_path_arca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_477917events_gtis.pkl"
    #gti_path_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_1787352events_gtis.pkl"

    # --- File paths as before ---

    if source == "vela":
        source_file = "/home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_analysis/inputs/sources/Vela_X-1.h5"
        Emin = 2e5
        Emax = 2e6
        source_name = "Vela"
        frequency = 11
        
        if angle == 30:
            input_file = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_478124events_vela30deg.hdf5"
            gti_path_arca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_478124events_gtis_vela30deg.pkl"
            input_file_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_1788128events_vela30deg.hdf5"
            gti_path_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_3090430events_gtis_crab30deg.pkl"
        elif angle == 25:
            input_file = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_328121events_vela25deg.hdf5"
            gti_path_arca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_328121events_gtis_vela25deg.pkl"
            input_file_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_1248520events_vela25deg.hdf5"
            gti_path_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_1248520events_gtis_vela25deg.pkl"
        elif angle == 20:
            input_file = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_203257events_vela20deg.hdf5"
            gti_path_arca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_203257events_gtis_vela20deg.pkl"
            input_file_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_845659events_vela20deg.hdf5"
            gti_path_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_845659events_gtis_vela20deg.pkl"
        elif angle == 15:
            input_file = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_112272events_vela15deg.hdf5"
            gti_path_arca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_112272events_gtis_vela15deg.pkl"
            input_file_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_507690events_vela15deg.hdf5"
            gti_path_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_507690events_gtis_vela15deg.pkl"
        elif angle == 10:
            input_file = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_49261events_vela10deg.hdf5"
            gti_path_arca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_49261events_gtis_vela10deg.pkl"
            input_file_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_232614events_vela10deg.hdf5"
            gti_path_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_232614events_gtis_vela10deg.pkl"
        elif angle == 5:
            input_file = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_12943events_vela5deg.hdf5"
            gti_path_arca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_12943events_gtis_vela5deg.pkl"
            input_file_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_58211events_vela5deg.hdf5"
            gti_path_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_58211events_gtis_vela5deg.pkl"
    elif source == "crab":
        source_file = "/home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_analysis/inputs/sources/CrabPulsar.h5"
        Emin = 7e4
        Emax = 6e5
        source_name = "Crab"
        frequency = 30
        if angle == 30:
            input_file = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_1958155events_crab30deg.hdf5"
            gti_path_arca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_1958155events_gtis_crab30deg.pkl"
            input_file_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_3090430events_crab30deg.hdf5"
            gti_path_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_3090430events_gtis_crab30deg.pkl"
        elif angle == 5:
            input_file = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_88478events_crab5deg.hdf5"
            gti_path_arca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_88478events_gtis_crab5deg.pkl"
            input_file_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_261929events_crab5deg.hdf5"
            gti_path_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_261929events_gtis_crab5deg.pkl"




    #source_file = "/home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_analysis/inputs/sources/CrabPulsar.h5"
    #source_name = "Crab"
    #input_file = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_12858events.hdf5"
    #input_file = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/arca_KM3NeT_00000133_0001_combined_477917events.hdf5"
    
    #input_file = "/home/wecapstor3/capn/capn107h/combined_eventlists/arca_KM3NeT_00000133_0001_combined_4982323events_mc.hdf5"
    #input_file_orca = "/home/wecapstor3/capn/capn107h/combined_eventlists/orca_KM3NeT_00000148_0001_combined_1008510events_mc.hdf5"
    #input_file_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_58211events.hdf5"
    #input_file_orca = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists/orca_KM3NeT_00000148_0001_combined_1787352events.hdf5"
    
    aeff_file = "/home/hpc/capn/capn107h/software/psr/src/scripts/aeff_summary_v9_flavors_arca.hdf5"
    aeff_file_orca = "/home/hpc/capn/capn107h/software/psr/src/scripts/aeff_summary_v9_flavors_orca.hdf5"

    shower_histo_file  = "/home/hpc/capn/capn107h/software/psr/src/scripts/arca_energy_response_histogram_arca_shower_pdgid12_-12_16_-16.hdf5"
    track_histo_file = "/home/hpc/capn/capn107h/software/psr/src/scripts/arca_energy_response_histogram_arca_track_pdgid14_-14.hdf5"
    shower_histo_file_orca  = "/home/hpc/capn/capn107h/software/psr/src/scripts/orca_energy_response_histogram_orca_shower_pdgid12_-12_16_-16.hdf5"
    track_histo_file_orca = "/home/hpc/capn/capn107h/software/psr/src/scripts/orca_energy_response_histogram_orca_track_pdgid14_-14.hdf5"

    file_path_angres_sh = "/home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_analysis/inputs/angular_res_files/AngularResolutionOverEnergy_aashower_mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013.hdf5"
    file_path_angres_tr = "/home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_analysis/inputs/angular_res_files/AngularResolutionOverEnergy_jmuon_mcv8.1.gsg_numu-CCHEDIS_1e2-1e8GeV.sirene.jterbr000132.hdf5"
    file_path_angres_sh_orca = "/home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_analysis/inputs/angular_res_files/AngularResolutionOverEnergy_jshower_gsg_elec-CC_1.0-100.0GeV_orca.hdf5"
    file_path_angres_tr_orca = "/home/hpc/capn/capn107h/software/psr/workflows/km3net/km3net_analysis/inputs/angular_res_files/AngularResolutionOverEnergy_jmuon_gsg_muon-CC_1.0-100.0GeV_orca.hdf5"

    if fluxtype == "classic":
        selected_detectors = ['arca'] #, 'orca', 'combined']
        E0 = 1e0
        phi0_values = np.logspace(-6, 5, 25)
    elif fluxtype == "bounded":
        selected_detectors = ['arca']
        E0 = Emin
        phi0_values = np.logspace(-14, -6, 25)

    #phi0_values = np.logspace(-14, -6, 25)
    #phi0_values = np.logspace(-6, 5, 5)

    gamma_values = np.linspace(1.5, 3, 25)
    

    # --- Load detector configurations ---
    detectors = {}
    gti_arrays = {}

    if 'arca' in selected_detectors or 'combined' in selected_detectors:
        detectors['arca'] = load_static_inputs(
            input_file, source_file, aeff_file, 
            shower_histo_file, track_histo_file, 
            file_path_angres_sh, file_path_angres_tr, 
            'arca'
        )
        gti_arrays['arca'] = np.array(loadGTIs(gti_path_arca)) if gti_path_arca else None

    if 'orca' in selected_detectors or 'combined' in selected_detectors:
        detectors['orca'] = load_static_inputs(
            input_file_orca, source_file, aeff_file_orca, 
            shower_histo_file_orca, track_histo_file_orca, 
            file_path_angres_sh_orca, file_path_angres_tr_orca, 
            'orca'
        )
        gti_arrays['orca'] = np.array(loadGTIs(gti_path_orca)) if gti_path_orca else None

    if 'combined' in selected_detectors:
        detectors['combined'] = None
        # Merge GTIs if both are available
        if gti_arrays.get('arca') is not None and gti_arrays.get('orca') is not None:
            gti_arrays['combined'] = np.vstack([gti_arrays['arca'], gti_arrays['orca']])
        elif gti_arrays.get('arca') is not None:
            gti_arrays['combined'] = gti_arrays['arca']
        elif gti_arrays.get('orca') is not None:
            gti_arrays['combined'] = gti_arrays['orca']
        else:
            gti_arrays['combined'] = None

    # --- Parameters ---
    #phi0_values = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]
    #phi0_values =  [3e4, 2e4, 1e4, 3e3, 2e3, 1e3, 3e2, 2e2, 1e2, 3e1, 2e1, 1e1, 3, 2, 1, 3e-1, 2e-1, 1e-1, 3e-2, 2e-2, 1e-2, 3e-3, 2e-3, 1e-3, 3e-4, 2e-4, 1e-4]
    #phi0_values =  [3e4, 2e4, 1e4, 3e3, 2e3, 1e3, 3e2, 2e2, 1e2, 3e1, 2e1, 1e1, 3, 2, 1, 3e-1, 2e-1, 1e-1, 3e-2, 2e-2, 1e-2, 3e-3, 2e-3, 1e-3, 3e-4, 2e-4, 1e-4,3e-5,2e-5, 1e-5, 3e-6, 2e-6, 1e-6]
    #gamma_values = [1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5 ] #, 2.0, 3.0]
    #gamma_values = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5]
    #phi0_values = [1e-12, 2e-12, 5e-12, 7e-12, 1e-11, 2e-11, 5e-11, 7e-11, 1e-10, 2e-10, 3e-10, 1e-9, 2e-9, 3e-9, 1e-8, 2e-8, 3e-8]
    #gamma_values = [1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0]
    #reference_starttime = 1661990400  # 2022-09-01 UTC

    #phi0_values = [8.1633e-9]
    gamma_values= [2.0]
    phi0_values = [1.4e-1]
    #length_values = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150]  # days


    #phi0_values = np.logspace(-14, -6, 25)
    #phi0_values = np.logspace(-6, 5, 5)
    #gamma_values = np.linspace(1.5, 3, 5)

    reference_starttime = 1.68e9
    results = []

    #E0 = 1e0
    #E0 = Emin

    #frequency = 11.0
    nbins = 32

    #max_angle = 30.0
    length_values  = [150.0]
    n_repeats = 1
    #inspect_combo = (5e-12,1.0)
    inspect_combo = (1e-1, 2.0)
    epoch_folding_distributions = {}
    hist_data = {}

    # Output directory for inspection plots
    plot_dir = f"{datetime_str}_inspection_plots_phi0_{inspect_combo[0]:.1e}_gamma_{inspect_combo[1]:.2f}"
    os.makedirs(plot_dir, exist_ok=True)

    total_jobs = len(phi0_values) * len(gamma_values)
    job_counter = 0

    total_runs = len(phi0_values) * len(gamma_values) * len(selected_detectors)* len(length_values) * n_repeats
    global_run_counter = 0

    

    for phi0 in phi0_values:
        for gamma in gamma_values:
            for length in length_values:
                job_counter += 1
                print(f"\n=== Job {job_counter}/{total_jobs} for phi0={phi0:.1e}, gamma={gamma:.2f} ===")

                for detector_label in selected_detectors:

                    print(f"-- Running detector: {detector_label.upper()} --")
                    max_ef_list = []
                    efstat_list_all_runs = []  # only stored for inspect_combo

                    for repeat_idx in range(n_repeats):
                        global_run_counter += 1
                        try:
                            if detector_label == 'combined':
                                cut_arca = run_detector_analysis(phi0, gamma, E0, detectors['arca'], frequency, angle, length,reference_starttime,Emin=Emin,Emax=Emax,source_name=source_name, zenith_cut = zenith_cut, fluxtype = fluxtype)
                                cut_orca = run_detector_analysis(phi0, gamma, E0, detectors['orca'], frequency, angle, length,reference_starttime,Emin=Emin,Emax=Emax,source_name=source_name, zenith_cut = zenith_cut, fluxtype = fluxtype)
                                combined = vstack([cut_arca, cut_orca])
                            else:
                                combined = run_detector_analysis(phi0, gamma, E0,  detectors[detector_label], frequency, angle, length,reference_starttime,Emin=Emin,Emax=Emax,source_name=source_name, zenith_cut = zenith_cut, fluxtype = fluxtype)

                            print(f"-----> Cuts Applied to {detector_label}...")
                            print(f"len(combined): {len(combined)}")

                            #print("combined: ", combined)
                            gti_array = gti_arrays.get(detector_label, None)
                            #print("gti_array: ", gti_array)

                            subset, used_live_time, used_gtis = extract_subset_within_gtis(combined, gti_array, length)
                            #print("subset: ",subset)
                            final_times = np.array(TimeSeries(time=Time(subset['time'], format='unix'))['time'].value)
                            final_times = np.sort(final_times)

                            if len(final_times) < 10:
                                max_ef = np.nan
                            else:
                                freqs_t = get_testfrequencies(frequency, 500, 2e-9)
                                freqs , efstat = epoch_folding_search(
                                    final_times,
                                    freqs_t,
                                    nbin=nbins,
                                    segment_size=1e10,
                                    gti=used_gtis
                                )
                                max_ef = np.max(efstat)

                            #print("-----> Eoch Folding Completed...")
                            #print("#################################")
                            #print("#                               #")
                            #print("Max Efstat: ", np.max(efstat))
                            #print(f"Phi0 {phi0} - gamma {gamma}")
                            #print(f"Detector {detector_label}")
                            #print(f"Combination {job_counter}/{total_jobs}")
                            #print(f"Iteration: {repeat_idx}/{n_repeats}")
                            #print("#                               #")
                            #print("#################################")

                            box_width = 50
                            title = " Epoch Folding Completed "
                            border = "#" * box_width
                            spacer = f"# {' ' * (box_width - 2)}#"

                            print(border)
                            print(f"# {title:^{box_width-4}} #")
                            print(border)
                            print(f"# Max EF stat : {np.max(efstat):<28}#")
                            print(f"# Phi0        : {phi0:<28}#")
                            print(f"# Gamma       : {gamma:<28}#")
                            print(f"# Detector    : {detector_label:<28}#")
                            print(f"# Length      : {length:<28}#")
                            print(f"# Combination : {job_counter}/{total_jobs:<18}#")
                            print(f"# Iteration   : {repeat_idx+1}/{n_repeats:<18}#")
                            print(f"# Overall run : {global_run_counter}/{total_runs:<18}#")
                            print(border)
                            print()

                            max_ef_list.append(max_ef)

                            # Save detailed distribution if it's the inspected combo
                            if (phi0, gamma) == inspect_combo:
                                efstat_list_all_runs.append((freqs,efstat))

                        except Exception as e:
                            print(f"    Failed for {detector_label}: {e}")
                            max_ef_list.append(np.nan)

                                    # Store plotting data
                if (phi0, gamma) == inspect_combo:
                    epoch_folding_distributions[detector_label] = efstat_list_all_runs
                    hist_data[detector_label] = max_ef_list
                    #print("Inspection Parameter Combination: EF Info attatched for later Plotting.")

                    #    results.append({'phi0': phi0, 'gamma': gamma, 'detector': detector_label, 'max_efstat': max_ef})
                    #except Exception as e:
                    #    print(f"    Failed for {detector_label}: {e}")
                    #    results.append({'phi0': phi0, 'gamma': gamma, 'detector': detector_label, 'max_efstat': np.nan})

                # Compute mean and error
                mean_max_ef = np.nanmean(max_ef_list)
                mean_error = np.nanstd(max_ef_list) / np.sqrt(np.sum(~np.isnan(max_ef_list)))

                print("#############################################################")
                print("-------------------------------------------------------------")
                print("                                                             ")
                print("                                                             ")
                print(f"               mean_max_ef: {mean_max_ef}                   ")
                print("                                                             ")
                print("                                                             ")
                print("-------------------------------------------------------------")
                print("#############################################################")
                
                results.append({
                    'zenith_cut': zenith_cut,
                    'angle': angle,
                    'phi0': phi0,
                    'gamma': gamma,
                    'detector': detector_label,
                    'length_days': length,
                    'mean_max_efstat': mean_max_ef,
                    'mean_max_efstat_error': mean_error
                })



            # Save all epoch folding distributions
            for det, all_runs in epoch_folding_distributions.items():
                plt.figure(figsize=(8,6))
                for freqs, ef in all_runs:
                    plt.plot(freqs, ef, alpha=0.3)
                plt.title(f"Epoch Folding Distributions - {det}")
                plt.xlabel(r"Frequency [Hz]")
                plt.ylabel(r"Epoch Folding $\chi^2$")

                plt.axhline(nbins - 1, ls='dotted', lw=1, color='k', label="Expected no signal")
                plt.axvline(frequency, ls='--', lw=1, alpha=0.5, color='k', label='Test frequency')
                plt.legend()

                plt.tight_layout()  
                plt.savefig(os.path.join(plot_dir, f"{det}_epoch_folding_distributions_smfig_bounded.png"), dpi=300)
                plt.close()

            # Save histogram of max_ef values
            for det, max_vals in hist_data.items():
                clean_vals = [v for v in max_vals if not np.isnan(v)]
                plt.figure(figsize=(8, 6))
                
                # Plot histogram
                plt.hist(clean_vals, bins=15, alpha=0.7, color='skyblue', edgecolor='black')
                plt.title(f"Maximum Chi2 Histogram - {det}")
                plt.xlabel(r"Maximum $\chi^2$")
                plt.ylabel("Count")

                # Calculate mean and error
                mean_val = np.nanmean(clean_vals)
                mean_err = np.nanstd(clean_vals) / np.sqrt(len(clean_vals))
                
                # Add mean line
                plt.axvline(mean_val, color='gold', linestyle='--', lw=2, label=f"Mean = {mean_val:.2f}")
                
                # Add shaded error band
                plt.axvspan(mean_val - mean_err, mean_val + mean_err, color='gold', alpha=0.2,
                            label=f"$\pm$ Error = {mean_err:.2f}")

                plt.legend()
                plt.tight_layout()
                plt.savefig(os.path.join(plot_dir, f"{det}_max_ef_histogram_smfig_bounded.png"), dpi=300)
                plt.close()

            #print(f"Inspection plots saved in: {plot_dir} , or maybe also not, who knows?")



    

    df = pd.DataFrame(results)
    #df.to_csv(f"sensitivity_study/final4/efstat_grid_results_withzenithcut{zenith_cut}_{fluxtype}flux_{source_name}_{angle}deg_{n_repeats}iterations.csv", index=False)
    df.to_csv(f"efstat_grid_results_withzenithcut{zenith_cut}_{fluxtype}flux_{source_name}_{angle}deg_{n_repeats}iterations_Emin1e4GeV.csv", index=False)
    print("Saved results.")


if __name__ == "__main__":
    main()