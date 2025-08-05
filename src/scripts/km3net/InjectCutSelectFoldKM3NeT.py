"""Inject an artificial signal and apply cuts and time selection and epoch folding algorithm for sensitivity study.
Usage: Sensitivity_Study_KM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR -e ENERGY_PDF [--rate=<float>] [--pulseshape=<pulseshape>] [--df=<float>] [--frequency=<float>] [--baseline=<float>] [--a=<float>] [--phi=<float>] [--kappa=<float>] [--plot=<plot>] [--method=<method>] [--E_min=<float>] [--E_max=<float>] [--gamma=<float>] [--E_cut=<float>] [--detector=<detector>] [--energy_min=<float>] [--energy_max=<float>] [--nhits_tr_min=<float>] [--nhits_sh_min=<float>] [--lik_tr_min=<float>] [--lik_sh_min=<float>] [--beta0_tr_max=<float>] [--zenith_min=<float>] [--trackscore_low_max=<float>] [--trackscore_high_min=<float>] [--gti_files=<gti_files>...] [--length=<float>] [--fraction=<fraction>] [--number_of_testf=<number_of_testf>] [--nbin=<nbin>] [--testdf=<float>] [--ratio=<ratio>] [--iteration=<iteration>] [--segment_size=<segment_size>] [--run_name=<run_name>]

Options:
  -h --help                                 Help
  -i --input_files INPUT_FILES              Input files
  -o --output_dir OUTPUT_DIR                Output file
  -e --energy_pdf ENERGY_PDF                Filepath of the hdf5 file storing the E_true - E_reco distribution
     --rate=<float>                         Rate of injected signal neutrinos in 1/d. [default: 1.]
     --pulseshape=<pulseshape>              Shape of the injected signal (sine or mvm). [default: mvm]
                                             if 'sine': 'df', 'frequency', 'baseline', 'a', 'phi' should be set
     --df=<float>                           Time resolution of the signal. [default: 0.1]
     --frequency=<float>                    Frequency of the signal. [default: 1.0]
     --baseline=<float>                     Offset on the y-axis. [default: 0.]
     --a=<float>                            Amplitude of the signal. [default: 1.]
     --phi=<float>                          Phase of the signal. [default: 0.]
     --kappa=<float>                        Shape parameter of the MVMD. [default: 5.]
     --plot=<plot>                          Bool, whether to plot the injected signal or not [default: False]
     --method=<method>                      String: 'classic' or 'base'. Method of Injection to use. [default: base]
     --E_min=<float>                        Minimum Energy of injected power law spectrum [default: 1e2]
     --E_max=<float>                        Max Energy of the spectrum [default: 1e8]
     --gamma=<float>                        Spectral index of the spectrum [default: 2.0]
     --E_cut=<float>                        Cutoff Energy of spectrum [default: 1e6]
     --detector=<string>                    Detector location ('arca','orca','antares') [default: arca]
     --energy_min=<float>                   Minimum energy in GeV [default: 0.]
     --energy_max=<float>                   Maximum energy in GeV [default: 1e10]
     --nhits_tr_min=<float>                 Minimum number of hits for track reconstructions [default: 0.]
     --nhits_sh_min=<float>                 Minimum number of hits for shower reconstructions [default: 0.]
     --lik_tr_min=<float>                   Minimum likelihood for track reconstructions [default: 0.]
     --lik_sh_min=<float>                   Minimum likelihood for shower reconstructions, negated [default: 0.]
     --beta0_tr_max=<float>                 Maximum beta0 value for track reconstructions [default: 1.0]
     --zenith_min=<float>                   Minimum zenith value in degrees (Below horizon: >90deg) [default: 0.]
     --trackscore_low_max=<float>           Maximum trackscore to keep at the lower edge of the score distribution [default: 0.5]
     --trackscore_high_min=<float>          Minimum trackscore to keep at the higher end [default: 0.5]
     --gti_files=<gti_files>...             Optional GTI files
     --length=<float>                       Length of the subsate to be created in days, if '--fraction' is set to True this has to be a fraction of the original length between 0 and 1 [default: 1.]
     --fraction=<fraction>                  Bool, whether a fractino of the time or an absolute time is set [default: False]
     --number_of_testf=<int>                Number of testfrequencies to test around the principle frequency. [default: 200]
     --testdf=<float>                       Resolution of testfrequencies. [default: 1e-5]
     --nbin=<int>                           Number of bins in the folded profile. [default: 32]
     --iteration=<int>                      Nr of current iteration [default: 0]
     --segment_size=<float>                 Length of the segments to be averaged in the periodogram [default: 5000]
     --run_name=<string>                    Name of the run [default: sensitivity]
""" 

from docopt import docopt
import os, glob
import h5py as h5py
import re
import numpy as np
from astropy.time import Time, TimeDelta
from astropy.timeseries import TimeSeries
import astropy.units as u
import stingray
from plens.PulseModel import MVMD, sinusoid

# PLENS Imports
import plens.EventList as EL
import plens.antares_hdf5
import plens.antares_hdf5 as antares_hdf5


from stingray.pulse.search import epoch_folding_search, z_n_search, search_best_peaks
from epochfolding.stingray_epochfolding import savehdf5, get_testfrequencies
from epochfolding.gtis import loadGTIs, saveGTIs

from stingray import EventList, Lightcurve
import warnings
import matplotlib.pyplot as plt
import traceback
import sys
import inspect

import km3io.definitions as kd

def create_binned_light_curve(times, bin_size_seconds=60):
    """Plot a binned light curve for the given event times."""
    if len(times) == 0:
        print("No events to plot.")
        return
    
    start_time, end_time = np.min(times), np.max(times) 
    num_bins = int((end_time - start_time) / bin_size_seconds)
    
    event_counts, bin_edges = np.histogram(times, bins=num_bins, range=(start_time, end_time))
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    
    return bin_centers, event_counts

def plot_simulated_signal(lc_original, new_event_times, frequency, output_file, num_cycles=5, bin_size_seconds=None):
    """
    Plots the pure simulated signal and overlays the actual injected event times.
    
    Parameters:
        lc_original : Lightcurve
            Light curve of the pure simulated signal.
        new_event_times : np.array
            Injected event times (overlaid as vertical lines).
        num_cycles : int, optional
            Number of signal periods to display (default is 5).
    """


    plt.figure(figsize=(10, 5))

    # Calculate period from frequency
    # frequency = 1 / np.median(np.diff(lc_original.time))  # Approximate frequency
    period = 1 / frequency  # Period of the signal
    if bin_size_seconds == None: bin_size_seconds = period/20.0

    # Define the time range to plot
    t_min = lc_original.time[0] + period
    t_max = t_min + (num_cycles + 2) * period

    # Mask to select data in the desired time range
    mask = (lc_original.time >= t_min) & (lc_original.time <= t_max)

    plt.plot(lc_original.time[mask], lc_original.counts[mask]*2, lw=3, label="Pulse Timing Model", color='darkorange', alpha = 0.7, zorder = 3)

    new_event_times_bin_centers, new_event_times_counts = create_binned_light_curve(new_event_times,bin_size_seconds=bin_size_seconds)
    plt.plot(new_event_times_bin_centers, new_event_times_counts, drawstyle='steps-mid', lw=2, label="Injected Events Lightcurve", color='firebrick', alpha=0.7, zorder = 3)

    # Overlay injected event times within the time range
    for event in new_event_times:
        if t_min <= event <= t_max:
            plt.axvline(event, color='darkgrey', alpha=0.6,  label='Injected Events' if event == new_event_times[0] else None, zorder = 2)

    

    

    plt.xlabel("Time (s)")
    plt.ylabel("Counts per Bin")
    plt.title(f"Simulated Signal with Injected Event Times ({num_cycles} cycles)")
    plt.legend(['Simulated Signal', 'Injected Events Lightcurve', 'Injected Events'])
    plt.xlim([t_min+0.5*period,t_max-1.5*period])
    plt.grid(alpha=0.5)

    output_plot = output_file + "simulated_signal.png"
    plt.savefig(output_plot)
    print(f"Plot saved: {output_plot}")


def plot_data_comparison(original_times, injected_signal_times, frequency, output_file, rate, num_cycles=5, bin_size_seconds = None):
    """
    Plots a comparison between the initial dataset and the injected light curve.
    
    Parameters:
        lc_initial : Lightcurve
            Light curve of the original data before injection.
        lc_injected : Lightcurve
            Light curve of the dataset after injection.
        frequency : float
            The frequency of the signal (used to calculate the period).
        num_cycles : int, optional
            Number of signal periods to display (default is 5).
    """
    plt.figure(figsize=(10, 5))

    
    
    # Calculate period from frequency
    period = 1 / frequency  # Period of the signal
    if bin_size_seconds == None: bin_size_seconds = period/20.0

    # Define the time range to plot
    t_min = np.min([original_times[0],injected_signal_times[0]]) + period
    t_max = t_min + (num_cycles + 2) * period

    

    original_times_bin_centers, original_event_counts = create_binned_light_curve(original_times,bin_size_seconds=bin_size_seconds)
    injected_signal_bin_centers, injected_signal_counts = create_binned_light_curve(injected_signal_times,bin_size_seconds=bin_size_seconds)

    

    # Plot both initial and injected light curves within the specified time range
    plt.plot(original_times_bin_centers, original_event_counts, drawstyle='steps-mid', lw=4, label="Original Data", color='darkorange', alpha=1)
    plt.plot(injected_signal_bin_centers, injected_signal_counts, drawstyle='steps-mid', lw=2, label="Injected Data", color='darkblue', alpha=1)
    
    plt.xlabel("Time (s)")
    plt.ylabel("Intensity")
    plt.title(f"Initial vs Injected Light Curve ({num_cycles} cycles - {frequency:.2e} Hz - Rate {rate} 1/d)")
    plt.xlim([t_min+0.5*period,t_max-1.5*period])
    plt.legend()
    plt.grid(alpha=0.5)

    output_plot = output_file + "signal_comparison.png"
    plt.savefig(output_plot)
    print(f"Plot saved: {output_plot}")

def sample_power_law_energies(n, E_min, E_max, gamma, E_cut):
    """
    Samples `n` energies from a power-law distribution: flux(E) ∝ E^(-gamma)
    between E_min and E_max.
    """
    if gamma == 1.0:
        # Special case: integral diverges, use logarithmic sampling
        r = np.random.uniform(0, 1, n)
        energies = E_min * (E_max / E_min) ** r
    else:
        # Inverse transform sampling
        r = np.random.uniform(0, 1, n)
        exponent = 1.0 - gamma
        E_min_pow = E_min ** exponent
        E_max_pow = E_max ** exponent
        energies = (E_min_pow + (E_max_pow - E_min_pow) * r) ** (1.0 / exponent)
    
    return energies

def sample_power_law_energies_cutoff(n, E_min, E_max, gamma, E_cut):
    """
    Vectorized sampling from a power-law with exponential cutoff.
    Returns exactly `n` samples between E_min and E_max.
    """
    energies = []
    max_pdf = (E_min ** -gamma) * np.exp(-E_min / E_cut)
    batch_size = max(1000, n)  # ensures efficiency for small n too

    while len(energies) < n:
        # Step 1: Draw trial samples from uniform proposal
        E_trial = np.random.uniform(E_min, E_max, size=batch_size)

        # Step 2: Compute target PDF values for these samples
        pdf_vals = (E_trial ** -gamma) * np.exp(-E_trial / E_cut)

        # Step 3: Acceptance probability is ratio to maximum
        accept_prob = pdf_vals / max_pdf

        # Step 4: Accept or reject
        accepted = E_trial[np.random.rand(batch_size) < accept_prob]
        energies.extend(accepted.tolist())

    return np.array(energies[:n])



def injectSignalRedistributeHighRes(time, energy, rate, n_events_inject, bin_time, pulseshape, frequency, baseline, a, phi, hist2d, bins_true, bins_reco, E_min=1.1e5, E_max=2e6, E_cut=1e6, gamma=2.0, kappa=None, plot=False):
    """
    Injects a signal with a given pulseshape into an event list.

    Parameters:
        time : np.array
            Original event times.
        energy : np.array
            Original event energies
        rate: float
            Rate by which signal neutrinos should roughly be injected
        n_events_injects : int
            Number of injected signal neutrino events.
        bin_time : float
            Time bin width.
        pulseshape : str
            Type of pulse shape ('sine' or 'mvm').
        frequency : float
            Frequency of the pulse train.
        baseline : float
            Y-axis offset.
        a : float
            Amplitude.
        phi : float
            Phase.
        kappa : float, optional
            Shape parameter for MVM.
        plot : bool
            Whether to plot results.

    Returns:
        If plot=True:
            tuple: (lc_original, lc_injected, new_event_times, combined_events)
        Else:
            np.array: combined event times.
    """
    # Create a high-resolution time grid independent of `time`

    bin_time = 1e-1

    high_res_time = np.arange(time.min(), time.max(), bin_time)  # 10x finer resolution
    #print("len(high_res_time): ", len(high_res_time))
    #print("len(time): ", len(time))

    # Generate smooth signal counts based on the chosen pulse shape
    if pulseshape == 'mvm':
        # Use the MVMD function for the MVM pulse shape
        counts = MVMD(high_res_time, frequency, phi, kappa, a, baseline=baseline)
    elif pulseshape == 'sine':
        # Use the sinusoid function for the sine wave pulse shape
        counts = sinusoid(high_res_time, frequency, baseline, a, phi)

    # Generate a smooth light curve
    lc_original = Lightcurve(high_res_time, counts, dt=bin_time / 10, skip_checks=True)

    # Simulate event times from the smooth light curve
    ev = EventList()
    ev.simulate_times(lc_original)
    new_event_times = ev.time

    #new_event_energies = np.full_like(new_event_times, fill_value=1000.0, dtype=np.float64)
    new_event_energies_true = sample_power_law_energies(
        n=len(new_event_times),
        E_min=E_min,       # Adjust as needed
        E_max=E_max,     # Adjust as needed
        gamma=gamma,       # Slope of the spectrum
        E_cut=E_cut      # Energy cutoff
    )

    new_event_energies_reco, failed = sample_reco_energy_array(new_event_energies_true, hist2d, bins_true, bins_reco) 
    print(f"Successfully sampled: {(~failed).sum()} / {len(new_event_energies_reco)}")
    # Determine how many events to keep from the original event list
    total_events = len(time)

    

    #num_new_events = int(round(estimated_neutrinos * ratio))  # Number of new events to inject
    num_new_events = n_events_inject
    num_original_events_to_keep = total_events - num_new_events  # Remaining original events to keep

    if num_original_events_to_keep < 0:
        raise ValueError("Number of original events to keep is negative. Ensure estimated_neutrinos is reasonable.")

    # Randomly sample the desired number of injected events
    if len(new_event_times) >= num_new_events:
        indices_to_keep_new = np.random.choice(len(new_event_times), num_new_events, replace=False)
        new_event_times = new_event_times[indices_to_keep_new]
    else:
        indices_to_keep_new = np.random.choice(len(new_event_times), num_new_events, replace=True)
        new_event_times = new_event_times[indices_to_keep_new]

    # Keep a subset of original events
    indices_to_keep_original = np.random.choice(len(time), num_original_events_to_keep, replace=False)
    remaining_original_times = time[indices_to_keep_original]
    remaining_original_energies = energy[indices_to_keep_original]


    # Combine the original and injected event times
    combined_times = np.concatenate((remaining_original_times, new_event_times))
    combined_energies = np.concatenate((remaining_original_energies, new_event_energies_reco))
    # Sort by time to keep everything aligned
    sort_indices = np.argsort(combined_times)
    combined_times = combined_times[sort_indices]
    combined_energies = combined_energies[sort_indices]

    # Generate light curve for the injected event list
    time_bins = np.arange(time.min(), time.max(), bin_time)
    print("Time bins: ", time_bins)
    print("Length of time bins: ", len(time_bins))
    counts_combined, _ = np.histogram(combined_times, bins=time_bins)
    time_centers = (time_bins[:-1] + time_bins[1:]) / 2
    lc_injected = Lightcurve(time_centers, counts_combined, dt=bin_time, skip_checks=True)

    # Return the results
    if plot:
        return lc_original, lc_injected, new_event_times, new_event_energies_reco, combined_times, combined_energies
    else:
        return combined_times, combined_energies
        
def generate_base_lightcurve(frequency, start_time, bin_time, pulseshape, baseline, a, phi, kappa=None, num_periods=20):
    """
    Generates a high-resolution light curve for a short base interval.

    Parameters:
        frequency : float
            Frequency of the pulse train.
        bin_time : float
            Time bin width.
        pulseshape : str
            Type of pulse shape ('sine' or 'mvm').
        baseline : float
            Y-axis offset.
        a : float
            Amplitude.
        phi : float
            Phase.
        kappa : float, optional
            Shape parameter for MVM.
        num_periods : int
            Number of periods to include in the base interval.

    Returns:
        tuple: (base_event_times, base_interval_length)
    """
    base_interval_length = num_periods / frequency  # Total duration of the base interval
    high_res_time = np.arange(start_time, start_time + base_interval_length, bin_time)  # High-resolution time grid


    dtype_size = np.dtype(np.float64).itemsize  # Size of one element in bytes (usually 8 bytes)
    num_elements = int(base_interval_length / bin_time)  # Estimated number of elements
    estimated_memory_MB = (num_elements * dtype_size) / (1024**2)  # Convert bytes to MB

    print(f"Estimated memory required: {estimated_memory_MB:.2f} MB")
    print(f"Memory required for high_res_time: {high_res_time.nbytes / (1024**2):.2f} MB")

    print("high_res_time: ", high_res_time[:5])

    # Generate smooth signal counts based on the chosen pulse shape
    if pulseshape == 'mvm':
        counts = MVMD(high_res_time, frequency, phi, kappa, a, baseline=baseline)
    elif pulseshape == 'sine':
        counts = sinusoid(high_res_time, frequency, baseline, a, phi)
    
    # Generate a smooth light curve
    lc_base = Lightcurve(high_res_time, counts, dt=bin_time, skip_checks=True)

    #print("lc_base:", lc_base)
    

    # Simulate event times from the base light curve
    ev = EventList()
    ev.simulate_times(lc_base)
    base_event_times = ev.time

    #print("base_event_times: ", base_event_times)
    print("len(base_event_times): ", len(base_event_times))


    return lc_base, base_event_times, base_interval_length


def injectSignalWithBaseInterval(time, energy, rate, n_events_inject, bin_time, pulseshape, 
                                 frequency, baseline, a, phi, hist2d, bins_true, bins_reco, kappa=None, num_periods=20,  E_min=1.1e5, E_max=2e6, E_cut=1e6, gamma=2.0, plot=False):
    """
    Injects a signal using a periodic base interval approach.

    Parameters:
        time : np.array
            Original event times.
        ratio : float
            Fraction of estimated neutrino events to inject.
        n_events_inject : int
            Number of signal events to inject.
        bin_time : float
            Time bin width.
        pulseshape : str
            Type of pulse shape ('sine' or 'mvm').
        frequency : float
            Frequency of the pulse train.
        baseline : float
            Y-axis offset.
        a : float
            Amplitude.
        phi : float
            Phase.
        kappa : float, optional
            Shape parameter for MVM.
        num_periods : int
            Number of periods to use in the base interval.
        plot : bool
            Whether to plot results.

    Returns:
        np.array: Combined event times.
    """
    start_time = np.min(time)
    print("start_time: ", start_time)
    print("end_time: ", np.max(time))
    
    # Generate base light curve and get event times
    lc_base, base_event_times, base_interval_length = generate_base_lightcurve(
        frequency, start_time, bin_time, pulseshape, baseline, a, phi, kappa, num_periods
    )

    #print("lc_base_loaded: ", lc_base)
    print("len(base_event_times_loaded): ", len(base_event_times))
    print("base_interval_length: ", base_interval_length)

    # Determine how many events to inject
    total_events = len(time)
    #num_new_events = int(round(int(estimated_neutrinos) * ratio))  
    num_new_events = n_events_inject
    num_original_events_to_keep = total_events - num_new_events  

    print("total_events: ", total_events)
    print("num_new_events: ", num_new_events)

    if num_original_events_to_keep < 0:
        raise ValueError("Number of original events to keep is negative. Check estimated_neutrinos.")

    # Randomly select event times from the base and spread across dataset
    total_time_span = time.max() - time.min()
    num_intervals = int(np.ceil(total_time_span / base_interval_length))

    print(f"total_time_span: {total_time_span}, Type: {type(total_time_span)}")
    print(f"num_intervals: {num_intervals}, Type: {type(num_intervals)}")
    print(f"total_time_span/num_intervals: {total_time_span/num_intervals}")

    # Sample random base events and replicate across time
    chosen_base_times = np.random.choice(base_event_times, num_new_events, replace=True)
    random_intervals = np.random.randint(0, num_intervals, size=num_new_events)
    random_shifts = random_intervals * base_interval_length
    new_event_times = chosen_base_times + random_shifts


    # Randomly sample the desired number of injected events
    if len(new_event_times) >= num_new_events:
        indices_to_keep_new = np.random.choice(len(new_event_times), num_new_events, replace=False)
        new_event_times = new_event_times[indices_to_keep_new]
    else:
        indices_to_keep_new = np.random.choice(len(new_event_times), num_new_events, replace=True)
        new_event_times = new_event_times[indices_to_keep_new]

    # Keep a subset of original events
    indices_to_keep_original = np.random.choice(len(time), num_original_events_to_keep, replace=False)
    remaining_original_times = time[indices_to_keep_original]
    remaining_original_energies = energy[indices_to_keep_original]
    print("len(new_event_times): ", len(new_event_times))
    new_event_energies_true = sample_power_law_energies(
        n=len(new_event_times),
        E_min=E_min,       # Adjust as needed
        E_max=E_max,     # Adjust as needed
        gamma=gamma,       # Slope of the spectrum
        E_cut=E_cut      # Energy cutoff
    )
    print("len(new_event_energies_true): ", len(new_event_energies_true))
    
    new_event_energies_reco, failed = sample_reco_energy_array(new_event_energies_true, hist2d, bins_true, bins_reco) 
    print(f"Successfully sampled: {(~failed).sum()} / {len(new_event_energies_reco)}")

    print(f"chosen_base_times: {chosen_base_times[:5]}, Length: {len(chosen_base_times)}, Type: {type(chosen_base_times)}")
    print(f"Min: {np.min(chosen_base_times)}")
    print(f"Max: {np.max(chosen_base_times)}")
    print(f"Min (relative): {np.min(chosen_base_times)-start_time}")
    print(f"Max (relative): {np.max(chosen_base_times)-start_time}")    
    print(f"random_intervals: {random_intervals[:5]}, Length: {len(random_intervals)}, Type: {type(random_intervals)}")
    print(f"Min: {np.min(random_intervals)}")
    print(f"Max: {np.max(random_intervals)}")
    print(f"random_shifts: {random_shifts[:5]}, Length: {len(random_shifts)}, Type: {type(random_shifts)}")
    print(f"new_event_times: {new_event_times[:5]}, Length: {len(new_event_times)}, Type: {type(new_event_times)}")
    print(f"new_event_energies_reco: {new_event_energies_reco[:5]}, Length: {len(new_event_energies_reco)}, Type: {type(new_event_energies_reco)}")
    # Check if the periodicity is actually there:
    
    # Ensure events are within valid range
    #new_event_times = new_event_times[(new_event_times >= time.min()) & (new_event_times <= time.max())]
    #print(f"After Range Selection: new_event_times: {new_event_times}, Length: {len(new_event_times)}, Type: {type(new_event_times)}")
    print("new_event_times Min: ", np.min(new_event_times))
    print("new_event_times Max: ", np.max(new_event_times))
    # Keep a subset of original events
    indices_to_keep_original = np.random.choice(len(time), num_original_events_to_keep, replace=False)
    remaining_original_times = time[indices_to_keep_original]
    remaining_original_energies = energy[indices_to_keep_original]

    print(f"indices_to_keep_original: {indices_to_keep_original}, Length: {len(indices_to_keep_original)}, Type: {type(indices_to_keep_original)}")
    print(f"remaining_original_events: {remaining_original_times}, Length: {len(remaining_original_times)}, Type: {type(remaining_original_times)}")

    # Combine original and injected event times
    combined_times = np.sort(np.concatenate((remaining_original_times, new_event_times)))
    combined_energies = np.concatenate((remaining_original_energies, new_event_energies_reco))

    print("Length of combined_event_times:", len(combined_times))
    print("Length of combined_energies:", len(combined_energies))
    # Sort by time to keep everything aligned
    sort_indices = np.argsort(combined_times)
    combined_times = combined_times[sort_indices]
    combined_energies = combined_energies[sort_indices]
    print(f"combined_events: {combined_times}, Length: {len(combined_times)}, Type: {type(combined_times)}")
    lc_injected = []

    # Return results
    if plot:
        return lc_base, lc_injected, new_event_times , combined_times, combined_energies
    else:
        return combined_times, combined_energies

def load_histogram(file_path):
    with h5py.File(file_path, 'r') as f:
        hist2d = f['hist2d'][()]
        bins_true = f['bins_true'][()]
        bins_reco = f['bins_reco'][()]
    return hist2d, bins_true, bins_reco

def sample_reco_energy_array(true_energies, hist2d, bins_true, bins_reco):
    """
    Samples reconstructed energies for an array of true energies based on the 2D histogram.

    Parameters:
    - true_energies (array-like): True energies to sample from
    - hist2d (2D array): Histogram [true_energy_bin, reco_energy_bin]
    - bins_true (1D array): Bin edges for true energy
    - bins_reco (1D array): Bin edges for reco energy

    Returns:
    - reco_samples (np.ndarray): Sampled reconstructed energies
    - failed_mask (np.ndarray): Boolean mask indicating which samples failed
    """
    true_energies = np.asarray(true_energies)
    reco_samples = np.full_like(true_energies, fill_value=np.nan, dtype=float)
    bin_centers_reco = 0.5 * (bins_reco[1:] + bins_reco[:-1])

    for i, E_true in enumerate(true_energies):
        true_bin_idx = np.digitize(E_true, bins_true) - 1

        # Check if within valid histogram range
        if 0 <= true_bin_idx < hist2d.shape[0]:
            slice_hist = hist2d[true_bin_idx, :]
            if np.sum(slice_hist) > 0:
                pdf = slice_hist / np.sum(slice_hist)
                reco_samples[i] = np.random.choice(bin_centers_reco, p=pdf)

    # Optional: mask of failed samples (outside valid range or empty bins)
    failed_mask = np.isnan(reco_samples)
    return reco_samples, failed_mask

def apply_selection_cuts(event_list, detector, energy_min, energy_max, nhits_tr_min,
                         nhits_sh_min, lik_tr_min, lik_sh_min, beta0_tr_max,
                         zenith_min, trackscore_low_max, trackscore_high_min):
    """Apply a set of selection cuts to a KM3NeT event list."""

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
    zenith_min_rad = zenith_min * np.pi / 180
    event_list = event_list[event_list['zenith'] >= zenith_min_rad]

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

    return event_list

def extract_subset_within_gtis(event_list, gti_array, desired_length_days):
    """Extract a time-limited subset from event_list, using GTIs to accumulate real live-time."""
    subset_event_indices = []
    accumulated_live_time = 0.0
    used_gtis = []

    event_times = event_list['time']

    if gti_array is not None:
        for start, stop in gti_array:
            gti_duration = stop - start

            if accumulated_live_time + gti_duration > desired_length_days:
                remaining_time = desired_length_days - accumulated_live_time
                new_stop = start + remaining_time
                in_partial_gti = (event_times.value >= start) & (event_times.value < new_stop)
                selected = np.where(in_partial_gti)[0]
                subset_event_indices.extend(selected.tolist())
                used_gtis.append([start, new_stop])
                accumulated_live_time += remaining_time
                break
            else:
                in_gti = (event_times.value >= start) & (event_times.value <= stop)
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

def main():
    arguments = docopt(__doc__)

    data = {key.replace("-", ""): arguments[key] for key in arguments}

    input_files = []
    for pattern in data['input_files']:
        input_files.extend(glob.glob(pattern))
    input_files.sort()

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    plot = str(data['plot']).lower() == "true"

    print(f"DEBUG: Plot option is {plot} (type: {type(plot)})")

    method = str(arguments['--method'])
    length = float(arguments['--length'])
    fraction = str(arguments['--fraction']).lower() == "true"

    energy_pdf_path = data['energy_pdf']
    hist2d, bins_true, bins_reco = load_histogram(energy_pdf_path)

    # Check how gti_files are present

    gti_files = data.get('gti_files', [])
    if gti_files is None:
        gti_files = []
    else:
        gti_files = [file for file in gti_files]

    gti_files.sort()

    # Validate GTI file count
    if len(gti_files) == 1 and len(input_files) > 1:
        print("Warning: Only one GTI file provided for multiple input files. It will be applied to all input files.")
    elif len(gti_files) not in [0, 1, len(input_files)]:
        print("Error: Number of GTI files must be either 0, 1, or equal to the number of input files.")
        return
    


    for idx, file in enumerate(input_files):

        current_gti = None

        if gti_files:
            gti_file = gti_files[0] if len(gti_files) == 1 else gti_files[idx]

            gti_table = loadGTIs(gti_file)
            gti_table = np.array(gti_table)
            print("gti_table: ",gti_table)

            gti_start = gti_table[:,0]
            gti_stop = gti_table[:,1]
            current_gti = np.array([gti_start, gti_stop]).T
            print(f"GTIs: {current_gti}")


        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]

        output_file = f"{data['output_dir']}{file_name}_{str(data['run_name'])}_{data['frequency']}Hz_signalwithcuts_{str(length)}days"

        with h5py.File(file) as input_file:
            EventList = EL.readEventList(input_file)
            times = EventList['time'].value.astype(float)
            zenith = EventList['zenith'].value.astype(float)
            energy = EventList['energy'].value.astype(float)
            rate = float(data['rate'])

            print(f"Max Time: {times.max()}")
            print(f"Min Time: {times.min()}")
            duration = times.max() - times.min()
            duration_days = duration/(24.0*3600.0)
            print(f"Duration: {duration} s - {duration_days} d")
            n_events = len(times)
            print(f"Number of Events: {n_events}")
            rate_max = 1.0*n_events/duration_days
            print(f"Max Rate: {rate_max} 1/d")

            if rate > rate_max:
                print(f"Rate ({rate}) exceeds maximum possible rate ({rate_max})! That doesn't work, signal Injection is skipped...")
                sys.exit(1)
            else:
                print(f"Rate ({rate}) is lower than maximum possible rate ({rate_max}), all fine!")
                n_events_inject = int(rate*duration_days)
                print(f"Number of injected events: {n_events_inject}")

            inject_kwargs = {
                'time': times,
                'energy': energy,
                'rate': float(data['rate']),
                'n_events_inject' : n_events_inject,
                'bin_time': float(data['df']),
                'pulseshape': data['pulseshape'],
                'frequency': float(data['frequency']),
                'baseline': float(data['baseline']),
                'a': float(data['a']),
                'phi': float(data['phi']),
                'hist2d' : hist2d,
                'bins_true' : bins_true,
                'bins_reco' : bins_reco,
                'E_min': float(data['E_min']),
                'E_max': float(data['E_max']),
                'E_cut': float(data.get('E_cut', 1e6)),     # fallback to default if not provided
                'gamma': float(data.get('gamma', 2.0)),     # fallback to default if not provided   
            }
            if data['pulseshape'] == 'mvm':
                inject_kwargs['kappa'] = float(data['kappa'])
                output_file = output_file + '_mvm_Emin' + str(data['E_min']) + 'g' + str(data['gamma'])
            else:
                #output_file += '_sine.hdf5'
                output_file = output_file + '_sine_Emin' + str(data['E_min']) + 'g' + str(data['gamma'])

            

            # Only use simulated times from within the gtis when implemented

            if plot == True:
                try:
                    if method == 'classic':
                        lc_original, lc_injected, new_event_times, new_event_energies, combined_event_times , combined_energies= injectSignalRedistributeHighRes(**inject_kwargs, plot=plot)
                    elif method == 'base':
                        lc_original, lc_injected, new_event_times, combined_event_times, combined_energies = injectSignalWithBaseInterval(**inject_kwargs)
                    else: print("--method must be 'classic' or 'base'!")
                    # Try plotting
                    plot_simulated_signal(lc_original, new_event_times, frequency=float(data['frequency']), output_file=output_file)
                    plot_data_comparison(times, combined_event_times, frequency=float(data['frequency']), output_file=output_file, rate=float(data['rate']))

                except Exception as e:
                    print(f"Warning: Plotting failed due to error: {e}")
                    traceback.print_exc()  # Print full error details

                    # Generate dummy PNG files
                    dummy_files = [f"{output_file}simulated_signal.png", f"{output_file}signal_comparison.png"]
                    for dummy_file in dummy_files:
                        print(f"Creating dummy plot file: {dummy_file}")
                        with open(dummy_file, "w") as f:
                            f.write("Dummy plot file - plotting failed.")
            else:
                if method == 'classic':
                    combined_event_times, combined_energies = injectSignalRedistributeHighRes(**inject_kwargs)
                elif method == 'base':
                    combined_event_times, combined_energies = injectSignalWithBaseInterval(**inject_kwargs)
                else: print("--method must be 'classic' or 'base'!")

            TimeEventList= TimeSeries(time=Time(combined_event_times, format='unix'))
            InjectedEventList = EventList.copy()
            InjectedEventList['time'] = TimeEventList['time']
            InjectedEventList['energy'] = combined_energies #* EventList['energy'].unit

            print(f"Signal Injected. Frequency: {float(data['frequency'])} Hz ")

            cut_params = {
                'detector': str(arguments['--detector']),
                'energy_min': float(arguments['--energy_min']),
                'energy_max': float(arguments['--energy_max']),
                'nhits_tr_min': float(arguments['--nhits_tr_min']),
                'nhits_sh_min': float(arguments['--nhits_sh_min']),
                'lik_tr_min': float(arguments['--lik_tr_min']),
                'lik_sh_min': float(arguments['--lik_sh_min']),
                'beta0_tr_max': float(arguments['--beta0_tr_max']),
                'zenith_min': float(arguments['--zenith_min']),
                'trackscore_low_max': float(arguments['--trackscore_low_max']),
                'trackscore_high_min': float(arguments['--trackscore_high_min']),
            }

            InjectedEventListWithCuts = apply_selection_cuts(InjectedEventList, **cut_params)
            print(f"All cuts applied.")


            ######## SUBSET SELECTION ####################################

            event_times = InjectedEventListWithCuts['time']
            print("event_times[:10]: ", event_times[:10])
            print("current_gti[:3]: ", current_gti[:3])

            
            # Apply GTI filtering mask as safety check
            if current_gti is not None:
                mask = np.zeros_like(event_times, dtype=bool)
                for start, stop in current_gti:
                    in_gti = (event_times.value >= start) & (event_times.value <= stop)
                    mask |= in_gti
                filtered_events = InjectedEventListWithCuts[mask]
                filtered_times = event_times[mask]
            else:
                filtered_events = InjectedEventListWithCuts
                filtered_times = event_times

            print("filtered_times[:10]: ", filtered_times[:10])

            # Safety check
            if len(filtered_times) == 0:
                print(f"No events after GTI filtering in file {file}")
                continue

            # Compute total active time
            if current_gti is not None:
                total_active_time = np.sum(current_gti[:, 1] - current_gti[:, 0])
            else:
                total_active_time = filtered_times[-1] - filtered_times[0] # This messes up, not calculating time in seconds
            
            total_active_time = event_times[-1].value - event_times[0].value
            # Determine target subset length
            if fraction:
                if not (0 <= length <= 1):
                    print(f"Invalid fraction value: {length}. Must be between 0 and 1.")
                    continue
                desired_length = total_active_time * length
            else:
                desired_length = length * 24.0*3600.0

            total_active_time_days = total_active_time / (3600.0*24.0)
            print(f"Total Active Time: {total_active_time} s - {total_active_time_days} days")
            print(f"Chosen Length: {desired_length} s - {desired_length/(24.0*3600.0)} d")
            print(f"Number of Events: {len(event_times)}")

            # Extract the GTI-respecting subset
            InjectedEventListWithCutsSubset, used_live_time, used_gtis = extract_subset_within_gtis(InjectedEventListWithCuts, current_gti, desired_length)

            used_live_time_days = used_live_time /(24.0*3600.0)
            print(f"New Amount of Events for {used_live_time_days} days ({used_live_time:.2f} s): {len(InjectedEventListWithCutsSubset['time'])}")

            InjectedEventListWithCutsSubset.meta['length'] = length



            for key, val in cut_params.items():
                InjectedEventListWithCuts.meta[key] = val

        output_file_table = output_file + "table.hdf5"
        if os.path.exists(output_file_table):
            print("File already existed. Deleting File...")
            os.remove(output_file_table)

        write_data = False

        if write_data is True:
            with h5py.File(output_file_table, 'w') as output:
                InjectedEventListWithCutsSubset.write(output, format='hdf5', overwrite=True, serialize_meta=True)
            print(f"File written to: {output_file_table}")

        n_total_events = len(InjectedEventListWithCutsSubset['time'])
        total_events_file = os.path.join(str(data['output_dir']),f"{n_total_events}_aftercuts_subset_total_events.txt" )
        with open(total_events_file, "w") as f:
            f.write(str(n_total_events) + "\n") 
            f.write(str(n_events_inject) + "\n")
            f.write(str(length) + "\n")

        ########### EPOCH FOLDING ##############################################

        FinalEventList= TimeSeries(time=Time(InjectedEventListWithCutsSubset['time'], format='unix'))

        print(f"Times: {np.array(FinalEventList['time'].value)}")
        frequencies = get_testfrequencies(float(data['frequency']), int(data['number_of_testf']), float(data['testdf']))

        print(f"First 10 Times: {[f'{x:.3f}' for x in np.array(FinalEventList['time'].value)[:10]]}")
        print(f"Last 10 Times: {[f'{x:.3f}' for x in np.array(FinalEventList['time'].value)[-10:]]}")
        print(f"Number of available events: {len(np.array(FinalEventList['time'].value))}")
        print(f"Test frequencies: {frequencies}")
        print(f"nbins: {int(data['nbin'])}")

        freq, efstat = epoch_folding_search(
            np.array(FinalEventList['time'].value),
            frequencies,
            nbin=int(data['nbin']),
            segment_size=float(data['segment_size']),
            gti=used_gtis
        ) 

        output_file_ef = output_file + f"_r{data['rate']}_I{data['iteration'].zfill(4)}_{data['frequency']}Hz_ef.hdf5"
        output_plot_ef = output_file + f"_r{data['rate']}_I{data['iteration'].zfill(4)}_{data['frequency']}Hz_efplot.png"

        if plot is True:
            plt.figure()
            plt.plot(freq, efstat, label='EF statistics', alpha=0.8)
            plt.axhline(int(data['nbin']) - 1, ls='--', lw=3, color='k', label='n - 1')
            plt.axvline(float(data['frequency']), lw=3, alpha=0.5, color='r', label='Correct frequency')
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('EF Statistics')
            _ = plt.legend()

            threshold = (np.max(efstat)*0.1+int(data['nbin']) - 1)
            best_x, best_y = search_best_peaks(freq,efstat,threshold)
            y_min, y_max = plt.ylim()
            offset = 0.04 * (y_max - y_min)

            '''
            for i, (x_value, y_value) in enumerate(zip(best_x, best_y)):
                label = 'peaks' if i == 0 else None  # Label only the first line
                plt.axvline(x_value, ls='dotted', lw=2, color='k', label=label)

                # Annotate the peak with its y-value slightly to the righ
                #plt.text(x_value + 0.03e-5, y_value, f"{y_value:.2f}", color='darkorange', fontsize=10)
                plt.text(x_value + 0.03e-5, y_value, f"{x_value:.4e}", color='darkorange', fontsize=10)
                plt.text(x_value + 0.03e-5, y_value-offset, f"{(float(data['frequency'])/x_value):.2f}", color='darkorange', fontsize=10)
            '''
            plt.savefig(output_plot_ef)

            print(f"Plot saved: {output_plot_ef}")

        print(f"Stingray version: {stingray.__version__}")
        print(f"epoch_folding_search is defined in: {inspect.getfile(epoch_folding_search)}")

        with h5py.File(output_file_ef, 'w') as out:
            savehdf5(freq, efstat, out)


if __name__ == "__main__":
    main()
