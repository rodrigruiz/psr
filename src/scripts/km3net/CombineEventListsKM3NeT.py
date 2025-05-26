""" Combines multiple eventlist objects generated from KM3NeT data to a single one.

Usage: CombineEventListsKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR -s SOURCE_SPECS_FILE [--zenith_threshold=<zenith_threshold>] [--delta_search_min=<delta_search_min>] [--filestype=<filestype>] [--detector=<detector>]

Options:
  -h --help                              Show this help message
  -i --input_files INPUT_FILES...        Input files
  -o --output_dir OUTPUT_DIR             Output directory  
  -s --source SOURCE_SPECS_FILE          Hdf5 file containing information about the source of interest (ra, dec, P_orb, ...) 
     --zenith_threshold=<float>          Zenith threshold for events below the horizon [default: 90.0]
     --delta_search_min=<float>          Minimal angular search cone size in degrees [default: 8]
     --filestype=<string>                Type of the files ('data' or 'mc') [default: data]
     --detector=<string>                 Detector location ('arca','orca','arca&orca') [default: arca]
"""

# python3 psr/src/scripts/CombineEventListsKM3NeT.py -i "/home/hpc/capn/capn107h/software/correctedeventlistTestOutput/*" -o combinedeventlistTestOutput/


from docopt import docopt
import os
import glob
import h5py
import numpy as np
from astropy.table import vstack, Table
import astropy.units as u
import plens.EventList as EL

from stingray import EventList, Lightcurve
import warnings
import matplotlib.pyplot as plt
from astropy.time import Time
from h5py import File
from epochfolding.gtis import findGTIs, saveGTIs

import datetime

def estimate_total_neutrino_events_old(time, zenith_input, radian=False):
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
    # Step 1: Find indices where zenith > 90° (neutrino-only periods)
    if radian == False:
        zenith = zenith_input* u.rad.to(u.deg)
    else: zenith = zenith_input*u.deg
    below_horizon_mask = zenith > 90
    above_horizon_mask = zenith <= 90

    print("Zenith values:")
    print(zenith)  # Show the zenith angles
    
    print("Mask for zenith > 90° (neutrino-only):")
    print(below_horizon_mask)  # Show the mask for neutrino-only
    
    print("Mask for zenith <= 90° (muon-only):")
    print(above_horizon_mask)  # Show the mask for muon-only

    neutrino_only_times = time[below_horizon_mask]
    muon_only_times = time[above_horizon_mask]

    # Step 2: Compute total time spent in each condition
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

    return N_neutrino_scaled

def estimate_total_neutrino_events(time, zenith_input, zenith_threshold=90.0, n_events_min = 10 , duration_min=0*u.h, title_prefix=None, output_file=None, radian=False, plot_intervals=True, shade_intervals=True, print_rates = False):
    """
    Estimates the total number of neutrino events in the dataset by:
    1. Identifying time intervals where at least 10 consecutive events have zenith > 90°.
    2. Computing the event rate in those periods.
    3. Calculating the mean event rate across all valid intervals.
    4. Using this mean rate to estimate the total number of neutrino events.
    5. Returning a masked dataset of neutrino-only events for debugging.

    Parameters:
        time (np.array): Array of event times (in seconds).
        zenith_input (np.array): Array of zenith angles.
        radian (bool): If True, zenith_input is in radians, otherwise in degrees.
    
    Returns:
        tuple: (Estimated total neutrino events, Masked dataset for debugging)
    """
    print("Initial time array:", time)
    print("Initial zenith array:", zenith_input)
    print("Amount of Events Total: ", len(time))
    # Convert zenith angles to degrees if necessary
    zenith = zenith_input * u.rad.to(u.deg) if not radian else zenith_input * u.rad
    print("Converted zenith angles:", zenith)
    
    # Identify indices where zenith > 90°
    #below_horizon_mask = zenith > 90
    print("Zenith Angles: ", zenith)
    if radian: below_horizon_mask = zenith > zenith_threshold*np.pi/180.0 
    else: below_horizon_mask = zenith > zenith_threshold
    print(f"Boolean mask for zenith > {zenith_threshold}°:", below_horizon_mask)
    print("Amont of Masked Events: ", len(zenith[below_horizon_mask]))
    print("Masked Zenith Angles: ", zenith[below_horizon_mask])
    # Find start and end indices of continuous segments with zenith > 90°
    indices = np.where(below_horizon_mask)[0]
    print(f"Indices of events with zenith > {zenith_threshold}°:", indices)
    
    if len(indices) == 0:
        raise ValueError("No events found with zenith > 80°.")
    
    intervals = []  # Stores (start_time, end_time) of valid intervals
    event_rates = []  # Stores event rate for each interval
    masked_event_times = []  # Store events for debugging
    masked_event_zeniths = []

    start_idx = indices[0]
    for i in range(1, len(indices)):
        if indices[i] != indices[i - 1] + 1:
            # Found a break in continuity, process previous segment

            segment_length = indices[i - 1] - start_idx + 1
            start_time = time[start_idx]
            end_time = time[indices[i - 1]]
            duration = end_time - start_time

            print("duration: ", duration)
            print("duration_min: ", duration_min)


            if  segment_length >= n_events_min and duration >= duration_min:

                event_count = segment_length - 1  # +1
                event_rate = 1.0*(event_count)/duration # Check what makes sense here
                event_rates.append(event_rate)
                intervals.append((start_time, end_time))
                masked_event_times.extend(time[start_idx:indices[i - 1] + 1].value)
                masked_event_zeniths.extend(zenith[start_idx:indices[i - 1] + 1].value)
                print(f"Interval found: Start {start_time}, End {end_time}, Duration {duration}, Event Rate {event_rate}")
            start_idx = indices[i]  # Start new segment
    
    # Process last segment
    segment_length = indices[-1] - start_idx + 1
    start_time = time[start_idx]
    end_time = time[indices[-1]]
    duration = end_time - start_time
    
    if segment_length >= n_events_min and duration >= duration_min:

        event_count = segment_length - 1  #+1
        event_rate = 1.0*(event_count)/ duration # Check what makes sense here
        event_rates.append(event_rate)
        intervals.append((start_time, end_time))
        masked_event_times.extend(time[start_idx:indices[-1] + 1].value)
        masked_event_zeniths.extend(zenith[start_idx:indices[-1] + 1].value)
        print(f"Last interval: Start {start_time}, End {end_time}, Duration {duration}, Event Rate {event_rate}")

    total_duration = time.max() - time.min()
    print("Total observation duration:", total_duration)
    
    if len(event_rates) == 0:
        print(f"Warning: No valid neutrino-only intervals with at least {n_events_min} consecutive events. Default value of 0.1*n_events is chosen.")
        estimated_total_neutrinos = 0.1*len(time)
        mean_event_rate = estimated_total_neutrinos/total_duration
    elif len(time) < np.mean(u.Quantity(event_rates)) * total_duration:
        estimated_total_neutrinos = 0.5*len(time)
        mean_event_rate = np.mean(u.Quantity(event_rates))
        print("Mean event rate:", mean_event_rate)
    else:
        # Compute mean event rate
        mean_event_rate = np.mean(u.Quantity(event_rates))
        print("Mean event rate:", mean_event_rate)
        print("Mean event rate:", mean_event_rate.to(u.Hz))
        estimated_total_neutrinos = mean_event_rate * total_duration

    print("Estimated total neutrino events:", estimated_total_neutrinos)

    masked_event_zeniths_deg = u.Quantity(masked_event_zeniths, u.deg) 

    #print("masked_event_zeniths_deg: ",masked_event_zeniths_deg)
    masked_event_times = Time(masked_event_times,format="unix")


    plt.figure(figsize=(10, 6))
    plt.scatter(np.array(time.mjd),np.array(zenith.value),s=2,label="All Events",color="darkblue")
    plt.scatter(np.array(masked_event_times.mjd),np.array(masked_event_zeniths_deg.to_value(u.deg)),s=2,label="Masked Events",color="darkorange")
    plt.axhline(90,ls="dashed",c='gray',alpha=0.7,label="Horizon")
    plt.axhline(zenith_threshold,c='k',ls="dotted",label="Zenith Threshold")

    if plot_intervals:
        for start, end in intervals:
            plt.axvline(start.mjd, linestyle="dotted", color="gray", alpha=0.8)
            plt.axvline(end.mjd, linestyle="dotted", color="gray", alpha=0.8)
    
    if print_rates:
        for (start, end), rate in zip(intervals, event_rates):
            mid_time = (start.mjd + end.mjd) / 2
            y_pos = zenith_threshold + 5  # Position text slightly above the threshold line
            plt.text(mid_time, y_pos, f"{rate*3600:.2e} events/h", rotation=90, fontsize=8,
                    verticalalignment='bottom', horizontalalignment='center', color='black', alpha=0.8)
    
    if shade_intervals:
        for start, end in intervals:
            plt.axvspan(start.mjd, end.mjd, color='orange', alpha=0.2)

    

    plt.xlabel("Time (MJD)")
    plt.ylabel("Zenith Angle (°)")
    plt.title(f"{title_prefix} Zenith Angle over Time \n {len(time)} Total Events - {int(estimated_total_neutrinos)} Estimated Neutrino Events (Zenith Threshold: {int(zenith_threshold)})")
    plt.grid()
    
    plt.gca().invert_yaxis()
    plt.legend()

    output_plot = output_file + "zenith_over_time.png"
    plt.savefig(output_plot)
    print(f"Plot saved: {output_plot}")
    

    
    return estimated_total_neutrinos, masked_event_times, masked_event_zeniths, mean_event_rate, intervals, event_rates

def plot_binned_light_curve_old(mjd_times, masked_event_times, mean_event_rate, estimated_neutrino_events, title_prefix, output_file, bin_size_seconds=60, plot_intervals=False, shade_intervals=False, highlight_intervals=True, zenith_threshold = 90.0):
    """Plot a binned light curve for the given event times."""
    if len(mjd_times) == 0:
        print("No events to plot.")
        return
    
    print("Mean Event Rate: ", mean_event_rate.value)
    print("Mean event Rate per Bin: ", mean_event_rate.to(u.Hz)*bin_size_seconds)

    
    start_time, end_time = np.min(mjd_times), np.max(mjd_times)
    bin_size_days = bin_size_seconds / (24 * 3600)  # Convert seconds to days
    num_bins = int((end_time - start_time) / bin_size_days)
    
    event_counts, bin_edges = np.histogram(mjd_times, bins=num_bins, range=(start_time, end_time))
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    
    plt.figure(figsize=(10, 6))
    plt.plot(bin_centers, event_counts, drawstyle='steps-mid', color='darkblue', label=f"Above Horizon Events (<{zenith_threshold}deg zenith)")

    if highlight_intervals:
        masked_event_counts, _ = np.histogram(masked_event_times.mjd, bins=bin_edges)
        plt.plot(bin_centers, masked_event_counts, drawstyle='steps-mid', color='darkorange', label=f"Below Horizon Events (<{zenith_threshold}deg zenith)")
    
    if plot_intervals:
        for start, end in zip(bin_edges[:-1], bin_edges[1:]):
            plt.axvline(start, linestyle="dotted", color="gray", alpha=0.8)
    
    if shade_intervals:
        for start, end in zip(bin_edges[:-1], bin_edges[1:]):
            plt.axvspan(start, end, color='orange', alpha=0.1)

    plt.axhline(mean_event_rate.to(u.Hz).value * bin_size_seconds, color="gray", linestyle="dashed", label="Mean Rate per Bin")
    
    plt.xlabel("Time (MJD)")
    plt.ylabel("Event Count")
    plt.title(f"{title_prefix} Binned Light Curve (bin size: {bin_size_seconds} s\n Zenith Threshold {zenith_threshold} ° - {len(mjd_times)} Total Events")
    plt.grid()
    plt.legend()
    
    output_plot = output_file + "lightcurve.png"
    plt.savefig(output_plot)
    print(f"Plot saved: {output_plot}")

def plot_binned_light_curve(mjd_times, masked_event_times, mean_event_rate, estimated_neutrino_events, title_prefix, output_file, bin_size_seconds=60, zenith_threshold=90.0):
    """Plot a binned light curve for the given event times with color-coded intervals."""
    if len(mjd_times) == 0:
        print("No events to plot.")
        return
    
    print("Mean Event Rate: ", mean_event_rate.value)
    print("Mean event Rate per Bin: ", mean_event_rate.to(u.Hz) * bin_size_seconds)

    start_time, end_time = np.min(mjd_times), np.max(mjd_times)
    bin_size_days = bin_size_seconds / (24 * 3600)  # Convert seconds to days
    num_bins = int((end_time - start_time) / bin_size_days)
    
    # Histogram for all events
    event_counts, bin_edges = np.histogram(mjd_times, bins=num_bins, range=(start_time, end_time))
    bin_centers = (bin_edges[1:] + bin_edges[:-1]) / 2
    
    # Histogram for below-horizon events
    masked_event_counts, _ = np.histogram(masked_event_times.mjd, bins=bin_edges)
    
    # Create a color array: orange for below-horizon intervals, blue otherwise
    colors = ["darkorange" if masked_event_counts[i] > 0 else "darkblue" for i in range(len(bin_centers))]
    
    plt.figure(figsize=(10, 6))
    
    # Plot the histogram using the color array
    for i in range(len(bin_centers) - 1):
        plt.plot([bin_centers[i], bin_centers[i + 1]], [event_counts[i], event_counts[i + 1]],
                 drawstyle='steps-mid', color=colors[i])
    
    # Mean rate per bin
    plt.axhline(mean_event_rate.to(u.Hz).value * bin_size_seconds, color="gray", linestyle="dashed", label="Mean Rate per Bin")
    
    plt.xlabel("Time (MJD)")
    plt.ylabel("Event Count")
    plt.title(f"{title_prefix} Binned Light Curve (bin size: {bin_size_seconds} s)\n Zenith Threshold {zenith_threshold}° - {len(mjd_times)} Total Events - {int(estimated_neutrino_events)} Neutrino Events")
    plt.grid()
    plt.legend()
    
    output_plot = output_file + "lightcurve.png"
    plt.savefig(output_plot)
    print(f"Plot saved: {output_plot}")

def event_energy_check(CombinedEventList, output_file):
    """
    Sorts the CombinedEventList by energy in descending order and saves the top 10 highest energy events to a file.
    
    Parameters:
        CombinedEventList (numpy structured array or similar format): The event list to be processed.
        output_file (str): The filename where the top 10 highest energy events will be saved.
    """
    # Sort by 'energy' in descending order
    CombinedEventList.sort("energy")
    CombinedEventList = CombinedEventList[::-1]  # Reverse to get highest energies first

    # Select the top 10 highest energy events
    top_10_events = CombinedEventList[:10]

    # Save to a text file
    with open(output_file, "w") as f:
        for event in top_10_events:
            f.write(str(event) + "\n")

    print(f"Top 10 highest energy events saved to {output_file}")

def event_rates_overview(intervals, event_rates, output_file):

    with open(output_file, 'w') as f:
        f.write("# Start_Time (mjd)\tEnd_Time (mjd)\tEvent_Rate\n")
        for (start, end), rate in zip(intervals, event_rates):
            f.write(f"{start}\t{end}\t{rate:.10f}\n")

def summarize_event_files(folder_path, file_suffix="total_events.txt", output_file="summary.txt"):
    summary_data = []

    for filename in os.listdir(folder_path):
        if filename.endswith(file_suffix):
            file_path = os.path.join(folder_path, filename)
            try:
                with open(file_path, 'r') as f:
                    lines = [line.strip() for line in f.readlines()]
                    if len(lines) < 4:
                        print(f"Skipping incomplete file: {filename}")
                        continue

                    n_total = int(lines[0])
                    n_estimated = float(lines[1])
                    dt_str = lines[2]
                    source_name = lines[3]
                    opening_angle = lines[4]
                    detector = lines[5]
                    filestype = lines[6]
                    

                    dt = datetime.datetime.strptime(dt_str, "%Y-%m-%d %H:%M:%S")
                    summary_data.append((dt, n_total, n_estimated, source_name, opening_angle, detector, filestype, filename))

            except Exception as e:
                print(f"Error processing {filename}: {e}")

    # Sort by datetime (descending: newest first)
    summary_data.sort(reverse=True, key=lambda x: x[0])

    # Write summary file
    #with open(os.path.join(folder_path, output_file), 'w') as out:
    #    out.write("# Datetime\t\t\tTotal_Events\tEstimated_Neutrinos\tSource\t\tFilename\n")
    #    for dt, total, estimated, source, fname in summary_data:
    #        out.write(f"{dt}\t\t{total}\t\t{estimated:.2f}\t\t\t{source}\t{fname}\n")

    with open(os.path.join(folder_path, output_file), 'w') as out:
        out.write(f"# {'Datetime':<20} {'Total_Events':>12} {'Estimated_Neutrinos':>20}  {'Source':<15} {'Min_Opening_Angle':<20} {'Detector':<20} {'File_Type':<15} Filename\n")
        for dt, total, estimated, source, opening_angle, detector, filestype, fname in summary_data:
            out.write(f"{dt.strftime('%Y-%m-%d %H:%M:%S')}  {total:>12}  {estimated:>20.2f}  {source:<15} {opening_angle:<20} {detector:<20} {filestype:<15} {fname}\n")

    print(f"Summary written to {output_file}")

def main():
    arguments = docopt(__doc__)

    input_files = arguments['--input_files']
    output_dir = arguments['--output_dir']
    zenith_threshold = float(arguments['--zenith_threshold'])
    detector = str(arguments['--detector'])

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with File(arguments['--source'], 'r') as f:
        source_name = f['source_name'][()].decode('utf-8')
        print("Source Name:", source_name)

    CombinedEventList = None
    time_intervals = []
    timediff_threshold = 1.0  

    for file in input_files:
        with h5py.File(file, 'r') as h5_file:
            EventList = EL.readEventList(h5_file)

            file_min = EventList['time'][0]
            file_max = EventList['time'][-1]
            time_intervals.append([file_min, file_max])

            if CombinedEventList is None:
                CombinedEventList = EventList
            else:
                CombinedEventList = vstack([CombinedEventList,EventList])
            
    time_intervals.sort(key=lambda x: x[0])

    # Merge intervals if they are close enough
    merged_intervals = []
    current_start, current_end = time_intervals[0]


    if time_intervals:
        current_start, current_end = time_intervals[0]

        for interval in time_intervals[1:]:
            next_start, next_end = interval

            if next_start - current_end < timediff_threshold:
                current_end = max(current_end, next_end)
            else:
                merged_intervals.append([current_start, current_end])
                current_start, current_end = next_start, next_end

        # Add the final interval
        merged_intervals.append([current_start, current_end])
    else:
        print("No time intervals found. Check input files.")
        merged_intervals = []

    # Add the last interval
    merged_intervals.append([current_start, current_end])
    print("Merged time intervals with uninterrupted data:")
    print(merged_intervals)
            
    CombinedEventList.sort("time")
    print("CombinedEventlist:")
    print(CombinedEventList)

    event_count = len(CombinedEventList)


    # Extract common prefix from input filenames
    common_prefix = os.path.commonprefix(input_files)
    # Remove any trailing non-alphanumeric characters from common prefix
    common_prefix = os.path.basename(common_prefix).rstrip("_-.")
    if not common_prefix:
        common_prefix = ""
    #output_file = os.path.join(output_dir, f"{common_prefix}_combined_eventlist.hdf5")
    output_file = os.path.join(output_dir, f"{detector}_{common_prefix}_combined_eventlist_{event_count}events.hdf5")
    gti_output_file = os.path.join(output_dir, f"{detector}_{common_prefix}_combined_eventlist_{event_count}events_gtis")
    
    saveGTIs(merged_intervals, gti_output_file)

    #times = CombinedEventList['time'] #.value.astype(float)
    zenith = CombinedEventList['zenith'] #.value.astype(float)

   

    event_times = CombinedEventList['time']
    astropy_event_times = Time(event_times, format="unix")
    mjd_times = astropy_event_times.mjd
    title_prefix = str(source_name) + " " + detector
    opening_angle_min = float(arguments['--delta_search_min'])
    filestype = str(arguments['--filestype'])
    

    
    estimated_neutrino_events, masked_times, masked_zeniths, mean_event_rate, intervals, event_rates = estimate_total_neutrino_events(astropy_event_times, zenith, zenith_threshold, n_events_min = 5, duration_min = 300*u.s, title_prefix=title_prefix, output_file=output_file)
    plot_binned_light_curve(mjd_times, masked_times, mean_event_rate, estimated_neutrino_events, title_prefix, output_file, bin_size_seconds=600,zenith_threshold=zenith_threshold)


    # Save the combined EventList
    CombinedEventList.write(output_file, format='hdf5', path='data', overwrite=True, serialize_meta = True)
    print(f"Combined EventList saved to {output_file}")

    n_total_events = len(event_times)
    total_events_file = os.path.join(output_dir,f"{n_total_events}_total_events.txt" )
    with open(total_events_file, "w") as f:
        f.write(str(n_total_events) + "\n") 
        f.write(str(estimated_neutrino_events) + "\n")
        f.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")
        f.write(str(source_name) + "\n")
        f.write(str(opening_angle_min) + "\n")
        f.write(str(detector) + "\n")
        f.write(str(filestype))

    energy_check_file = os.path.join(output_dir,f"{n_total_events}_energy_check.txt" )
    event_energy_check(CombinedEventList, energy_check_file)

    event_rates_overview_file = os.path.join(output_dir, f"{n_total_events}_event_rates_overview.txt")
    event_rates_overview(intervals, event_rates, event_rates_overview_file)

    #folder_path = "/home/hpc/capn/capn107h/software/nextflow_output/combined_eventlists"
    #summarize_event_files(folder_path)


    
if __name__ == "__main__":
    main()
    