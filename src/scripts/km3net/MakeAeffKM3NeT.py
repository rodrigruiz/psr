import numpy as np
import h5py
from astropy.table import Table
import numpy as np
import pandas as pd
import numba as nb
from numba import jit, njit, prange, float64, int64
from km3pipe.math import azimuth, zenith
import astropy.coordinates as ac
from astropy.time import Time
import astropy.units as u


def aeff_2D_eventlist(e_bins, t_bins, eventlist, gamma, n_runs, method=1):
    """
    Calculate the effective area in energy and zenith angle bins.

    Parameters
    ----------
    e_bins : Array
        of energy bins in GeV
    t_bins : Array
        of zenith angle bins in rad
    dataset : pandas.DataFrame
        with the events
    gamma : float, default 1.4
        spectral index of simulated events
    nevents : float, default 2e7
        number of generated events

    Returns
    -------
    2D-Array
        with the effective area in m^2 binned in energy and zenith angle bins

    """

    #if mcversion >= 9.0: T = 1.0
    #else: T = 365.25 * 24 *3600
        
    cos_theta = np.cos(eventlist["theta_mc"])
    theta_bins = np.searchsorted(t_bins, cos_theta) - 1

    #theta_bins = np.searchsorted(t_bins, eventlist['theta_mc']) - 1
    energy_bins = np.searchsorted(e_bins, eventlist['energy_mc']) - 1

    w2 = np.array(eventlist['normalized_weight'].value)
    E = np.array(eventlist['energy_mc'].value)
    
    if method == 1: aeff = fill_aeff_2D(e_bins, t_bins, energy_bins, theta_bins, w2, E, gamma)/n_runs
    elif method == 2: aeff = fill_aeff_2D_method2(e_bins, t_bins, energy_bins, theta_bins, w2, E)/n_runs
    else:
        raise ValueError("Unknown method. Choose method=1 or method=2.")
    return aeff


@njit(fastmath=False, parallel=True)
def fill_aeff_2D(e_bins, t_bins, energy_bins, theta_bins, w2, E, gamma):
    """
    numba accelerated helper function to calculate effective area.

    """
    #T = 365.25 * 24 * 3600
    #T = 1.0
    aeff = np.empty((len(e_bins) - 1, len(t_bins) - 1))
    for k in prange(len(e_bins) - 1):
        for i in range(len(t_bins) - 1):
            mask = (energy_bins == k) & (theta_bins == i)
            #d_omega = -(np.cos(t_bins[i + 1]) - np.cos(t_bins[i]))
            d_omega = (t_bins[i + 1] - t_bins[i])
            d_E = (e_bins[k + 1]) ** (1 - gamma) - (e_bins[k]) ** (1 - gamma)
            aeff[k, i] = (
                (1 - gamma)
                * np.sum(E[mask] ** (-gamma) * w2[mask])
                / (d_omega * d_E * 2 * np.pi)
            )

    return aeff


@njit(fastmath=False, parallel=True)
def fill_aeff_2D_method2(e_bins, t_bins, energy_bins, theta_bins, w2, E):
    """
    numba accelerated helper function to calculate effective area.

    """
    
    aeff = np.empty((len(e_bins) - 1, len(t_bins) - 1))
    for k in prange(len(e_bins) - 1):
        for i in range(len(t_bins) - 1):
            mask = (energy_bins == k) & (theta_bins == i)
            
            d_E = np.log10(e_bins[k + 1]/e_bins[k])
            d_Omega = 4 * np.pi/len(t_bins)
            #d_Omega = 2 * np.pi * (t_bins[i + 1] - t_bins[i])
            aeff[k, i] = (
                np.sum(E[mask] ** (-1) * w2[mask])
                / (np.log(10) * d_E * d_Omega)
            )

    return aeff


print("Functions built.")

def create_aeff_file_from_eventlist(config, output_path, E_min, E_max, gamma=1.4, method=1):
    """
    Computes and saves the 2D and 1D effective area for one or more neutrino flavors.

    Parameters
    ----------
    config : dict
        Dictionary with flavor-specific settings. Example:
        {
            "input_file": "/path/to/input.hdf5",
            "flavors": {
                "numu": {
                    "pdg_id": 14,
                    "nruns": 518,
                },
                "anue": {
                    "pdg_id": -12,
                    "nruns": 518,
                }
            }
        }
    output_path : str
        Path to save the resulting HDF5 file.
    gamma : float
        Spectral index of the simulated spectrum.
    method : int
        Method for computing the effective area (1 or 2).
    """

    #from your_module import aeff_2D_eventlist  # Replace with correct import

    print(f"Loading events from: {config['input_file']}")
    with h5py.File(config["input_file"], 'r') as h5_file:
        event_table = Table.read(h5_file)

    # Define binning
    #E_min, E_max = 1e2, 1e8
    n_energy_bins, n_theta_bins = 100, 100
    energy_bins = np.logspace(np.log10(E_min), np.log10(E_max), n_energy_bins + 1)
    cos_theta_bins = np.linspace(-1.0, 1.0, n_theta_bins + 1)
    dcos = np.diff(cos_theta_bins)

    # Output structure
    with h5py.File(output_path, "w") as f:
        f.create_dataset("energy_bins", data=energy_bins)
        f.create_dataset("cos_theta_bins", data=cos_theta_bins)
        f.attrs["gamma"] = gamma
        f.attrs["method"] = method

        for flavor, info in config["flavors"].items():
            print(f"Processing flavor: {flavor}")

            # Select events
            pdg_id = info["pdg_id"]
            nruns = info["nruns"]
            flavor_events = event_table[event_table["pdg_id"] == pdg_id]

            # Compute effective areas
            aeff_2d = aeff_2D_eventlist(
                energy_bins,
                cos_theta_bins,
                flavor_events,
                gamma=gamma,
                n_runs=nruns,
                method=method
            )

            aeff_1d = np.sum(aeff_2d * dcos[np.newaxis, :] * 0.5, axis=1)

            # Store in HDF5
            grp = f.create_group(flavor)
            grp.create_dataset("aeff_2d", data=aeff_2d)
            grp.create_dataset("aeff_1d", data=aeff_1d)
            grp.attrs["pdg_id"] = pdg_id
            grp.attrs["nruns"] = nruns

    print(f"Saved A_eff data to {output_path}")

def main():
    config_arca = {
        "input_file": "/home/wecapstor3/capn/capn107h/combined_eventlists/arca_KM3NeT_00000133_0001_combined_4983084events_mc.hdf5",
        "flavors": {
            "numu": {
                "pdg_id": 14,
                "nruns": 518
            },
            "anumu": {
                "pdg_id": -14,
                "nruns": 518
            },
            "anue": {
                "pdg_id": -12,
                "nruns": 518
            },
            "nue": {
                "pdg_id": 12,
                "nruns": 518
            },
            "nutau": {
                "pdg_id": 16,
                "nruns": 518
            },
            "anutau": {
                "pdg_id": -16,
                "nruns": 518
            }
        }
    }

    output_path_arca = "aeff_summary_v9_flavors_arca.hdf5"
    create_aeff_file_from_eventlist(config_arca, output_path_arca, 1e2, 1e8, gamma=1.4)

    config_orca = {
        "input_file": "/home/wecapstor3/capn/capn107h/combined_eventlists/orca_KM3NeT_00000148_0001_combined_1008510events_mc.hdf5",
        "flavors": {
            "numu": {
                "pdg_id": 14,
                "nruns": 481
            },
            "anumu": {
                "pdg_id": -14,
                "nruns": 481
            },
            "anue": {
                "pdg_id": -12,
                "nruns": 481
            },
            "nue": {
                "pdg_id": 12,
                "nruns": 481
            },
            "nutau": {
                "pdg_id": 16,
                "nruns": 481
            },
            "anutau": {
                "pdg_id": -16,
                "nruns": 481
            }
        }
    }

    output_path_orca = "aeff_summary_v9_flavors_orca.hdf5"
    create_aeff_file_from_eventlist(config_orca, output_path_orca, 1e0, 1e4, gamma=1.4)

if __name__ == "__main__":
    main()