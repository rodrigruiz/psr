"""Inject an artificial signal.
Usage: InjectSignal.py -i INPUT_FILES... -o OUTPUT_DIR [--ratio=<float>] [--pulseshape=<pulseshape>] [--df=<float>] [--frequency=<float>] [--baseline=<float>] [--a=<float>] [--phi=<float>] [--kappa=<float>] [--plot=<plot>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output file
     --ratio=<float>                     Ratio of injected signal to original count number. [default: 0.3]
     --pulseshape=<pulseshape>           Shape of the injected signal (sine or mvm). [default: mvm]
                                         if 'sine': 'df', 'frequency', 'baseline', 'a', 'phi' should be set
     --df=<float>                        Time resolution of the signal. [default: 0.1]
     --frequency=<float>                 Frequency of the signal. [default: 1]
     --baseline=<float>                  Offset on the y-axis. [default: 0.]
     --a=<float>                         Amplitude of the signal. [default: 1.]
     --phi=<float>                       Phase of the signal. [default: 0.]
     --kappa=<float>                     Shape parameter of the MVMD. [default: 5.]
     --plot=<plot>                       Bool, whether to plot the injected signal or not [default: False]
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


def injectSignalRedistributeHighRes(time, ratio, bin_time, pulseshape, frequency, baseline, a, phi, kappa=None, plot=False):
    """
    Injects a signal with a given pulseshape into an event list.

    Parameters:
        time : np.array
            Original event times.
        ratio : float
            Fraction of injected signal.
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
    high_res_time = np.arange(time.min(), time.max(), bin_time / 10)  # 10x finer resolution

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
    num_new_events = int(round(ratio * total_events))
    num_original_events_to_keep = total_events - num_new_events

    if num_original_events_to_keep < 0:
        raise ValueError("Number of original events to keep is negative. Ensure that ratio is between 0 and 1.")

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
    counts_combined, _ = np.histogram(combined_events, bins=time_bins)
    time_centers = (time_bins[:-1] + time_bins[1:]) / 2
    lc_injected = Lightcurve(time_centers, counts_combined, dt=bin_time, skip_checks=True)

    if plot:
        return lc_original, lc_injected, new_event_times, combined_events
    else:
        return combined_events


def main():
    arguments = docopt(__doc__)

    data = {key.replace("-", ""): arguments[key] for key in arguments}

    input_files = []
    for pattern in data['input_files']:
        input_files.extend(glob.glob(pattern))
    input_files.sort()

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    for file in input_files:
        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]

        output_file = f"{data['output_dir']}{file_name}_{data['frequency']}Hz_{data['a']}-signal"

        with h5py.File(file) as input_file:
            EventList = EL.readEventList(input_file)
            times = EventList['time'].value.astype(float)

            inject_args = [times, float(data['ratio']), float(data['df']), data['pulseshape'],
                           float(data['frequency']), float(data['baseline']), float(data['a']), float(data['phi'])]

            if data['pulseshape'] == 'mvm':
                inject_args.append(float(data['kappa']))
                output_file += '_mvm.hdf5'
            else:
                output_file += '_sine.hdf5'

            if data.get('plot', False):
                lc_original, lc_injected, new_event_times, combined_event_times = injectSignalRedistributeHighRes(*inject_args, plot=True)

                plot_simulated_signal(lc_original, new_event_times, frequency=float(data['frequency']) , output_file=output_file )
                plot_data_comparison(times, combined_event_times, frequency=float(data['frequency']) , output_file=output_file, ratio=float(data['ratio']))
            else:
                combined_event_times = injectSignalRedistributeHighRes(*inject_args)

            EventListNew = TimeSeries(time=Time(combined_event_times, format='unix'))
            InjectedEventList = EventList.copy()
            InjectedEventList['time'] = EventListNew['time']

        if os.path.exists(output_file):
            print("File already existed. Deleting File...")
            os.remove(output_file)

        with h5py.File(output_file, 'w') as output:
            InjectedEventList.write(output, format='hdf5', overwrite=True, serialize_meta=True)


if __name__ == "__main__":
    main()
