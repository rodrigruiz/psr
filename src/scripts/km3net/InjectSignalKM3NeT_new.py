"""Inject an artificial signal.
Usage: InjectSignal.py -i INPUT_FILES... -o OUTPUT_DIR --total_events_file TOTAL_EVENTS_FILE [--ratio=<float>] [--pulseshape=<pulseshape>] [--df=<float>] [--frequency=<float>] [--baseline=<float>] [--a=<float>] [--phi=<float>] [--kappa=<float>] [--plot=<plot>] [--method=<method>]

Options:
  -h --help                                 Help
  -i --input_files INPUT_FILES              Input files
  -o --output_dir OUTPUT_DIR                Output file
     --total_events_file TOTAL_EVENTS_FILE  File storing amount of total events used for timing analysis
     --ratio=<float>                        Ratio of injected signal to original count number. [default: 0.3]
     --pulseshape=<pulseshape>              Shape of the injected signal (sine or mvm). [default: mvm]
                                             if 'sine': 'df', 'frequency', 'baseline', 'a', 'phi' should be set
     --df=<float>                           Time resolution of the signal. [default: 0.1]
     --frequency=<float>                    Frequency of the signal. [default: 1]
     --baseline=<float>                     Offset on the y-axis. [default: 0.]
     --a=<float>                            Amplitude of the signal. [default: 1.]
     --phi=<float>                          Phase of the signal. [default: 0.]
     --kappa=<float>                        Shape parameter of the MVMD. [default: 5.]
     --plot=<plot>                          Bool, whether to plot the injected signal or not [default: False]
     --method=<method>                      String: 'classic' or 'base'. Method of Injection to use. [default: base]
""" 

from docopt import docopt
import os, glob
import h5py as h5py
import re
import numpy as np
from astropy.time import Time, TimeDelta
from astropy.timeseries import TimeSeries
import astropy.units as u
from plens.PulseModel import MVMD, sinusoid

# PLENS Imports
import plens.EventList as EL
import plens.antares_hdf5
import plens.antares_hdf5 as antares_hdf5

from stingray import EventList, Lightcurve
import warnings
import matplotlib.pyplot as plt
import traceback

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


def plot_data_comparison(original_times, injected_signal_times, frequency, output_file, ratio, num_cycles=5, bin_size_seconds = None):
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
    plt.title(f"Initial vs Injected Light Curve ({num_cycles} cycles - {frequency:.2e} Hz - SNR {ratio})")
    plt.xlim([t_min+0.5*period,t_max-1.5*period])
    plt.legend()
    plt.grid(alpha=0.5)

    output_plot = output_file + "signal_comparison.png"
    plt.savefig(output_plot)
    print(f"Plot saved: {output_plot}")


def injectSignalRedistributeHighRes_old(time, ratio, estimated_neutrinos, bin_time, pulseshape, frequency, baseline, a, phi, kappa=None, plot=False):
    """
    Injects a signal with a given pulseshape into an event list.

    Parameters:
        time : np.array
            Original event times.
        estimated_neutrinos : float
            Estimated number of neutrino events from previous calculations.
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
    high_res_time = np.arange(time.min(), time.max(), bin_time)  # 10x finer resolution
    print("len(high_res_time): ", len(high_res_time))
    print("len(time): ", len(time))

    # Generate smooth signal counts
    if pulseshape == 'mvm':
        counts = MVMD(high_res_time, frequency, phi, kappa, a, baseline=baseline)
    elif pulseshape == 'sine':
        counts = sinusoid(high_res_time, frequency, baseline, a, phi)

    # Generate a smooth light curve
    lc_original = Lightcurve(high_res_time, counts, dt=bin_time / 10, skip_checks=True)

    # Simulate event times from the smooth light curve
    ev = EventList()
    ev.simulate_times(lc_original)
    new_event_times = ev.time

    # Determine how many events to keep
    total_events = len(time)
    num_new_events = int(round(estimated_neutrinos * ratio))  # Use estimated neutrino count directly
    num_original_events_to_keep = total_events - num_new_events

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
    remaining_original_events = time[indices_to_keep_original]

    # Combine original and injected event times
    combined_events = np.sort(np.concatenate((remaining_original_events, new_event_times)))

    # Generate light curve for injected event list
    time_bins = np.arange(time.min(), time.max(), bin_time)
    print(time_bins)
    print(len(time_bins))
    counts_combined, _ = np.histogram(combined_events, bins=time_bins)
    time_centers = (time_bins[:-1] + time_bins[1:]) / 2
    lc_injected = Lightcurve(time_centers, counts_combined, dt=bin_time, skip_checks=True)

    if plot:
        return lc_original, lc_injected, new_event_times, combined_events
    else:
        return combined_events
    
def sample_power_law_energies(n, E_min, E_max, gamma, E_cut):
    """
    Sample energies from a power-law with exponential cutoff.

    Parameters:
        n : int
            Number of energies to sample.
        E_min : float
            Minimum energy.
        E_max : float
            Maximum energy.
        gamma : float
            Power-law index.
        E_cut : float
            Exponential cutoff energy.

    Returns:
        np.array of sampled energies.
    """
    from scipy.stats import rv_continuous
    from scipy.integrate import quad

    class PowerLawCutoff(rv_continuous):
        def _pdf(self, E):
            return E**(-gamma) * np.exp(-E / E_cut)

    # Normalize PDF
    dist = PowerLawCutoff(a=E_min, b=E_max, name='powerlaw_cutoff')
    norm = quad(dist._pdf, E_min, E_max)[0]
    dist._pdf = lambda E: (E**(-gamma) * np.exp(-E / E_cut)) / norm

    return dist.rvs(size=n)

    


def injectSignalRedistributeHighRes(time, energy, ratio, estimated_neutrinos, bin_time, pulseshape, frequency, baseline, a, phi, kappa=None, plot=False):
    """
    Injects a signal with a given pulseshape into an event list.

    Parameters:
        time : np.array
            Original event times.
        energy : np.array
            Original event energies
        estimated_neutrinos : float
            Estimated number of neutrino events from previous calculations.
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
    new_event_energies = sample_power_law_energies(
        n=len(new_event_times),
        E_min=1.1e5,       # Adjust as needed
        E_max=2e6,     # Adjust as needed
        gamma=2.0,       # Slope of the spectrum
        E_cut=1e6      # Energy cutoff
    )

    # Determine how many events to keep from the original event list
    total_events = len(time)
    num_new_events = int(round(estimated_neutrinos * ratio))  # Number of new events to inject
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
    combined_energies = np.concatenate((remaining_original_energies, new_event_energies))
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
        return lc_original, lc_injected, new_event_times, new_event_energies, combined_times, combined_energies
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


def injectSignalWithBaseInterval(time, ratio, estimated_neutrinos, bin_time, pulseshape, 
                                 frequency, baseline, a, phi, kappa=None, num_periods=20, plot=False):
    """
    Injects a signal using a periodic base interval approach.

    Parameters:
        time : np.array
            Original event times.
        ratio : float
            Fraction of estimated neutrino events to inject.
        estimated_neutrinos : float
            Estimated number of neutrino events from previous calculations.
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
    num_new_events = int(round(int(estimated_neutrinos) * ratio))  
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

    # Check if the periodicity is actually there:
    
    # Ensure events are within valid range
    #new_event_times = new_event_times[(new_event_times >= time.min()) & (new_event_times <= time.max())]
    #print(f"After Range Selection: new_event_times: {new_event_times}, Length: {len(new_event_times)}, Type: {type(new_event_times)}")
    print("new_event_times Min: ", np.min(new_event_times))
    print("new_event_times Max: ", np.max(new_event_times))
    # Keep a subset of original events
    indices_to_keep_original = np.random.choice(len(time), num_original_events_to_keep, replace=False)
    remaining_original_events = time[indices_to_keep_original]

    print(f"indices_to_keep_original: {indices_to_keep_original}, Length: {len(indices_to_keep_original)}, Type: {type(indices_to_keep_original)}")
    print(f"remaining_original_events: {remaining_original_events}, Length: {len(remaining_original_events)}, Type: {type(remaining_original_events)}")

    # Combine original and injected event times
    combined_events = np.sort(np.concatenate((remaining_original_events, new_event_times)))

    print(f"combined_events: {combined_events}, Length: {len(combined_events)}, Type: {type(combined_events)}")
    lc_injected = []

    # Return results
    if plot:
        return lc_base, lc_injected, new_event_times , combined_events
    else:
        return combined_events
    
def estimate_total_neutrino_events(time, zenith_input, radian=False):
    """
    Estimates the total number of neutrino events in the dataset by:
    1. Identifying time intervals where zenith > 90° (neutrino-only periods).
    2. Computing the event rate in those periods.
    3. Scaling that rate to the entire observation time.
    4. Computing the event rate for zenith ≤ 90° (atmospheric muon periods).
    
    Parameters:
        time (np.array): Array of event times (in seconds).
        zenith (np.array): Array of zenith angles (in degrees).
    
    Returns:
        float: Estimated total number of neutrino events in the dataset.
    """
    # Find indices where zenith > 90° (neutrino-only periods)
    if radian == False:
        zenith = zenith_input* u.rad.to(u.deg)
    else: zenith = zenith_input*u.deg
    below_horizon_mask = zenith > 90
    above_horizon_mask = zenith <= 90

    #print("Zenith values:")
    #print(zenith) 
    
    #print("Mask for zenith > 90° (neutrino-only):")
    #print(below_horizon_mask) 
    
    #print("Mask for zenith <= 90° (muon-only):")
    #print(above_horizon_mask)

    neutrino_only_times = time[below_horizon_mask]
    muon_only_times = time[above_horizon_mask]

    #Compute total time spent in each condition
    if len(neutrino_only_times) > 1:
        T_neutrino = np.sum(np.diff(neutrino_only_times))  # Sum of time gaps
    else:
        raise ValueError("Not enough neutrino-only events detected.")
    
    if len(muon_only_times) > 1:
        T_muon = np.sum(np.diff(muon_only_times))  # Sum of time gaps
    else:
        raise ValueError("Not enough muon-only events detected.")

    # Step 3: Compute total observation time
    T_total = time.max() - time.min()

    N_total = len(time)
    N_neutrino = len(neutrino_only_times)
    R_neutrino = N_neutrino / T_neutrino   # Neutrino event rate (events per second)

    N_muon = len(muon_only_times)
    R_muon = N_muon / T_muon  # Muon event rate (events per second)

    # Step 5: Scale up to the full observation time
    N_neutrino_scaled = R_neutrino * T_total

    print(f"Total neutrino-only observation time: {T_neutrino:.2f} s")
    print(f"Total atmospheric muon observation time: {T_muon:.2f} s")
    print(f"Total dataset observation time: {T_total:.2f} s")
    print(f"Neutrino event rate: {R_neutrino:.6f} events/s")
    print(f"Muon event rate: {R_muon:.6f} events/s")
    print(f"Estimated total neutrino events (scaled): {N_neutrino_scaled:.0f}")
    print(f"Total events: {N_total}")
    print(f"Ratio neutrino/all: {N_neutrino_scaled/N_total:.5f}")

    # Count events where zenith > 90° (neutrino-only) and zenith ≤ 90° (contains atmospheric muons)
    num_above_90 = np.sum(zenith > 90)
    num_below_90 = np.sum(zenith <= 90)
    # Compute the ratio
    ratio_above_below = num_above_90 / num_below_90 if num_below_90 > 0 else np.inf

    # Print the results
    print(f"Number of events with zenith > 90°: {num_above_90}")
    print(f"Number of events with zenith ≤ 90°: {num_below_90}")
    print(f"Ratio (above 90° / below 90°): {ratio_above_below:.4f}")

    return N_neutrino_scaled, neutrino_only_times



def main():
    arguments = docopt(__doc__)

    data = {key.replace("-", ""): arguments[key] for key in arguments}

    input_files = []
    for pattern in data['input_files']:
        input_files.extend(glob.glob(pattern))
    input_files.sort()

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    plot = data['plot'].lower() == "true"

    print(f"DEBUG: Plot option is {plot} (type: {type(plot)})")

        # Load total number of events from the file
    try:
        with open(arguments['--total_events_file'], "r") as f:
            lines = f.readlines()  # Read all lines into a list
            total_events = int(lines[0].strip())  # First line: total events (integer)
            estimated_neutrino_events = int(float(lines[1].strip()))  # Second line: estimated neutrino events (float)
    except Exception as e:
        print(f"Error reading total events file: {e}")
        total_events = None  # Default to None if reading fails
        estimated_neutrino_events = None  # Default to None if reading fails

    print(f"Total events used for timing analysis: {total_events}")
    print(f"Estimated neutrino events: {estimated_neutrino_events}")
    method = str(arguments['--method'])

    for file in input_files:
        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]

        output_file = f"{data['output_dir']}{file_name}_{data['frequency']}Hz_{data['a']}-signal"

        with h5py.File(file) as input_file:
            EventList = EL.readEventList(input_file)
            times = EventList['time'].value.astype(float)
            zenith = EventList['zenith'].value.astype(float)
            energy = EventList['energy'].value.astype(float)

            #estimated_neutrino_events = estimate_total_neutrino_events(times, zenith)

            #inject_args = [times, float(data['ratio']), float(data['df']), data['pulseshape'],
            #               float(data['frequency']), float(data['baseline']), float(data['a']), float(data['phi'])]
            
            inject_args = [times, energy, float(data['ratio']), estimated_neutrino_events, float(data['df']), data['pulseshape'],
               float(data['frequency']), float(data['baseline']), float(data['a']), float(data['phi'])]

            if data['pulseshape'] == 'mvm':
                inject_args.append(float(data['kappa']))
                output_file += '_mvm.hdf5'
            else:
                output_file += '_sine.hdf5'

            # Only use simulated times from within the gtis when implemented

            if plot == True:
                try:
                    if method == 'classic':
                        lc_original, lc_injected, new_event_times, new_event_energies, combined_event_times , combined_energies= injectSignalRedistributeHighRes(*inject_args, plot=plot)
                    elif method == 'base':
                        lc_original, lc_injected, new_event_times, combined_event_times = injectSignalWithBaseInterval(*inject_args)
                    else: print("--method must be 'classic' or 'base'!")
                    # Try plotting
                    plot_simulated_signal(lc_original, new_event_times, frequency=float(data['frequency']), output_file=output_file)
                    plot_data_comparison(times, combined_event_times, frequency=float(data['frequency']), output_file=output_file, ratio=float(data['ratio']))

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
                    combined_event_times, combined_energies = injectSignalRedistributeHighRes(*inject_args)
                elif method == 'base':
                    combined_event_times = injectSignalWithBaseInterval(*inject_args)
                else: print("--method must be 'classic' or 'base'!")

            TimeEventList= TimeSeries(time=Time(combined_event_times, format='unix'))
            InjectedEventList = EventList.copy()
            InjectedEventList['time'] = TimeEventList['time']
            InjectedEventList['energy'] = combined_energies #* EventList['energy'].unit

        if os.path.exists(output_file):
            print("File already existed. Deleting File...")
            os.remove(output_file)

        with h5py.File(output_file, 'w') as output:
            InjectedEventList.write(output, format='hdf5', overwrite=True, serialize_meta=True)

    n_total_events = len(times)
    total_events_file = os.path.join(str(data['output_dir']),f"{n_total_events}_total_events2.txt" )
    with open(total_events_file, "w") as f:
        f.write(str(n_total_events) + "\n") 
        f.write(str(estimated_neutrino_events) + "\n")


if __name__ == "__main__":
    main()
