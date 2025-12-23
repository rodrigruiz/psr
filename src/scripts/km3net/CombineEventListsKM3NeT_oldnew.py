""" Combines multiple eventlist objects generated from KM3NeT data to a single one.

Usage: CombineEventListsKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR -s SOURCE_SPECS_FILE [--zenith_threshold=<zenith_threshold>] [--delta_search_min=<delta_search_min>] [--filestype=<filestype>] [--detector=<detector>] [--combinedet=<combinedet>]

Options:
  -h --help                              Show this help message
  -i --input_files INPUT_FILES...        Input files
  -o --output_dir OUTPUT_DIR             Output directory  
  -s --source SOURCE_SPECS_FILE          Hdf5 file containing information about the source of interest (ra, dec, P_orb, ...) 
     --zenith_threshold=<float>          Zenith threshold for events below the horizon [default: 90.0]
     --delta_search_min=<float>          Minimal angular search cone size in degrees [default: 8]
     --filestype=<string>                Type of the files ('data' or 'mc') [default: data]
     --detector=<string>                 Detector location ('arca','orca','arca&orca') [default: arca]
     --combinedet=<string>               Whether this script is run for the combination of two detectors [default: False]
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

from astropy.time import Time
from epochfolding.gtis import saveGTIs

def compute_merged_intervals(input_files, timediff_threshold=1.0):
    """
    Read each event file and determine the earliest and latest event time
    (in UNIX seconds) from its event list. Then merge intervals if they
    are close enough.
    """
    time_intervals = []

    for file in input_files:
        with h5py.File(file, 'r') as h5_file:
            EventList = EL.readEventList(h5_file)

            if len(EventList) == 0:
                continue  # skip empty files

            # Min/max from raw event list (unfiltered)
            file_min = EventList['time'][0]
            file_max = EventList['time'][-1]

            time_intervals.append([file_min, file_max])

    if not time_intervals:
        print("No time intervals found. Check input files.")
        return []

    # Sort intervals by start time
    time_intervals.sort(key=lambda x: x[0])

    # Merge intervals if gaps are small
    merged_intervals = []
    current_start, current_end = time_intervals[0]

    for next_start, next_end in time_intervals[1:]:
        if next_start - current_end < timediff_threshold:
            current_end = max(current_end, next_end)
        else:
            merged_intervals.append([current_start, current_end])
            current_start, current_end = next_start, next_end

    # Append the last one
    merged_intervals.append([current_start, current_end])

    return merged_intervals


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
                    dt_str = lines[2]
                    source_name = lines[3]
                    opening_angle = lines[4]
                    detector = lines[5]
                    filestype = lines[6]
                    

                    dt = datetime.datetime.strptime(dt_str, "%Y-%m-%d %H:%M:%S")
                    summary_data.append((dt, n_total, source_name, opening_angle, detector, filestype, filename))

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
        out.write(f"# {'Datetime':<20} {'Total_Events':>12} {'Source':<15} {'Min_Opening_Angle':<20} {'Detector':<20} {'File_Type':<15} Filename\n")
        for dt, total, source, opening_angle, detector, filestype, fname in summary_data:
            out.write(f"{dt.strftime('%Y-%m-%d %H:%M:%S')}  {total:>12}  {source:<15} {opening_angle:<20} {detector:<20} {filestype:<15} {fname}\n")

    print(f"Summary written to {output_file}")

def main():
    arguments = docopt(__doc__)

    input_files = arguments['--input_files']
    output_dir = arguments['--output_dir']
    zenith_threshold = float(arguments['--zenith_threshold'])
    detector = str(arguments['--detector'])
    combinedet = str(arguments['--combinedet']).lower() == "true"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with File(arguments['--source'], 'r') as f:
        source_name = f['source_name'][()].decode('utf-8')
        print("Source Name:", source_name)

    CombinedEventList = None

    timediff_threshold = 1.0  

    for file in input_files:
        with h5py.File(file, 'r') as h5_file:
            EventList = EL.readEventList(h5_file)

            if CombinedEventList is None:
                CombinedEventList = EventList
            else:
                CombinedEventList = vstack([CombinedEventList,EventList])
            
            
    CombinedEventList.sort("time")
    print("CombinedEventlist:")
    print(CombinedEventList)

    event_count = len(CombinedEventList)


    # Extract common prefix from input filenames
    if combinedet:
        common_prefix = os.path.commonprefix(input_files)
        # Remove any trailing non-alphanumeric characters from common prefix
        common_prefix = os.path.basename(common_prefix).rstrip("_-.")
        if not common_prefix:
            common_prefix = ""
        #output_file = os.path.join(output_dir, f"{common_prefix}_combined_eventlist.hdf5")
        output_file = os.path.join(output_dir, f"{detector}_{common_prefix}_combinedet.hdf5")
        gti_output_file = os.path.join(output_dir, f"{detector}_{common_prefix}_combinedet_gtis")
    else:
        common_prefix = os.path.commonprefix(input_files)
        # Remove any trailing non-alphanumeric characters from common prefix
        common_prefix = os.path.basename(common_prefix).rstrip("_-.")
        if not common_prefix:
            common_prefix = ""
        #output_file = os.path.join(output_dir, f"{common_prefix}_combined_eventlist.hdf5")
        output_file = os.path.join(output_dir, f"{detector}_{common_prefix}_combined_{event_count}events.hdf5")
        gti_output_file = os.path.join(output_dir, f"{detector}_{common_prefix}_combined_{event_count}events_gtis")
    
    event_times = CombinedEventList['time']
    astropy_event_times = Time(event_times, format="unix")
    mjd_times = astropy_event_times.mjd
    title_prefix = str(source_name) + " " + detector
    opening_angle_min = float(arguments['--delta_search_min'])
    filestype = str(arguments['--filestype'])
    

    
    #estimated_neutrino_events, masked_times, masked_zeniths, mean_event_rate, intervals, event_rates = estimate_total_neutrino_events(astropy_event_times, zenith, zenith_threshold, n_events_min = 5, duration_min = 300*u.s, title_prefix=title_prefix, output_file=output_file)
    
    #plot_binned_light_curve(mjd_times, masked_times, mean_event_rate, estimated_neutrino_events, title_prefix, output_file, bin_size_seconds=600,zenith_threshold=zenith_threshold)


    # Save the combined EventList
    CombinedEventList.write(output_file, format='hdf5', path='data', overwrite=True, serialize_meta = True)
    print(f"Combined EventList saved to {output_file}")

    n_total_events = len(event_times)
    total_events_file = os.path.join(output_dir,f"{n_total_events}_total_events.txt" )
    with open(total_events_file, "w") as f:
        f.write(str(n_total_events) + "\n")
        f.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S") + "\n")
        f.write(str(source_name) + "\n")
        f.write(str(opening_angle_min) + "\n")
        f.write(str(detector) + "\n")
        f.write(str(filestype))

    
if __name__ == "__main__":
    main()
    