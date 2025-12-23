""" Select a subset of a certain time length from an EventList.  

Usage: SelectDataSetLenghtKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR [--gti_files=<gti_files>...] [--length=<float>] [--fraction=<fraction>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files or file pattern
  -o --output_dir OUTPUT_DIR             Output directory
     --gti_files=<gti_files>...          Optional GTI files
     --length=<float>                    Length of the subsate to be created in days, if '--fraction' is set to True this has to be a fraction of the original length between 0 and 1 
     --fraction=<fraction>               Bool, whether a fractino of the time or an absolute time is set [default: False]
"""
#python3 psr/src/scripts/CreateEventListKM3NeT.py -i '/home/hpc/capn/capn107h/software/hdf5TestOutput/*' -o eventlistTestOutput/ -s hdf5SourceFiles/Vela_X-1.h5

from docopt import docopt
import os, glob

from h5py import File, Group
from astropy.table import Table
from astropy.time import Time
from astropy.io.misc.hdf5 import read_table_hdf5, write_table_hdf5
import astropy.units as u
from astropy.coordinates import SkyCoord
import km3io.definitions as kd
import plens.EventList as EL
from km3astro.io import load_hdf5_tables

import numpy as np

from km3astro.coord import local_event
from km3astro import sources
from epochfolding.gtis import loadGTIs, saveGTIs

def extract_subset_within_gtis(event_list, event_times, gti_array, desired_length_days):
    """Extract a time-limited subset from event_list, using GTIs to accumulate real live-time."""
    subset_event_indices = []
    accumulated_live_time = 0.0
    used_gtis = []

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

def extract_subset_within_gtis_old(event_list, event_times, gti_array, desired_length_days):
    """Extract a time-limited subset from event_list, using GTIs to accumulate real live-time."""

    subset_event_indices = []
    accumulated_live_time = 0.0

    if gti_array is not None:
        # Iterate through GTIs
        for start, stop in gti_array:
            gti_duration = stop - start

            # Check how much time we still need
            if accumulated_live_time + gti_duration > desired_length_days:
                remaining_time = desired_length_days - accumulated_live_time
                new_stop = start + remaining_time
                in_partial_gti = (event_times.value >= start) & (event_times.value < new_stop)
                selected = np.where(in_partial_gti)[0]
                subset_event_indices.extend(selected.tolist())
                accumulated_live_time += remaining_time
                break
            else:
                # Full GTI fits
                in_gti = (event_times.value >= start) & (event_times.value <= stop)
                selected = np.where(in_gti)[0]
                subset_event_indices.extend(selected.tolist())
                accumulated_live_time += gti_duration

            if accumulated_live_time >= desired_length_days:
                break
    else:
        # No GTIs, just use a simple time cut from the start
        first_time = event_times[0]
        last_time = first_time + desired_length_days
        in_range = (event_times >= first_time) & (event_times < last_time)
        subset_event_indices = np.where(in_range)[0].tolist()
        accumulated_live_time = desired_length_days

    subset = event_list[subset_event_indices]
    return subset, accumulated_live_time

def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]

    input_files = []
    for pattern in data['input_files']:
        input_files.extend(glob.glob(pattern))
    input_files.sort()

    if not input_files:
        print(f"No files matching pattern: {input_files}")
        return
    
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

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    input_files = arguments['--input_files']
    output_dir = arguments['--output_dir']

    length = float(arguments['--length'])
    fraction = str(arguments['--fraction']).lower() == "true"

    for idx, file in enumerate(input_files):
        # Determine GTI for the current input file
        current_gti = None

        if gti_files:
            expocorr = True
            gti_file = gti_files[0] if len(gti_files) == 1 else gti_files[idx]

            #gti_table = read_table_hdf5(gti_file)
            #print(gti_file)
            gti_table = loadGTIs(gti_file)
            gti_table = np.array(gti_table)
            print("gti_table: ",gti_table)
            #print("Test")
            #print(gti_table[0])
            #print(gti_table[0][0].value)
            gti_start = gti_table[:,0]
            gti_stop = gti_table[:,1]
            current_gti = np.array([gti_start, gti_stop]).T
            print(f"GTIs: {current_gti}")

        with File(file, 'r') as h5_file:

            # Create Filename
            folder_path, file_name = os.path.split(file)
            file_name = os.path.splitext(file_name)[0]
            output_file = f"{arguments['--output_dir']}{file_name}_subset{str(length)}days.hdf5"


            # Load File:            
            EventList = EL.readEventList(h5_file)

            event_times = EventList['time']
            print("event_times[:10]: ", event_times[:10])
            
            # Apply GTI filtering mask (already done above when assigning 'filtered_events')
            if current_gti is not None:
                mask = np.zeros_like(event_times, dtype=bool)
                for start, stop in current_gti:
                    in_gti = (event_times.value >= start) & (event_times.value <= stop)
                    mask |= in_gti
                filtered_events = EventList[mask]
                filtered_times = event_times[mask]
            else:
                filtered_events = EventList
                filtered_times = event_times

            # Safety check
            if len(filtered_times) == 0:
                print(f"No events after GTI filtering in file {file}")
                continue

            # Compute total active time
            if current_gti is not None:
                print(current_gti[:,1])
                print(current_gti[:,0])
                print(current_gti[:,1]-current_gti[:,0])
                total_active_time = np.sum(current_gti[:, 1] - current_gti[:, 0])
            else:
                total_active_time = filtered_times[-1] - filtered_times[0]

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
            subset, used_live_time, used_gtis = extract_subset_within_gtis(filtered_events, filtered_times, current_gti, desired_length)

            used_live_time_days = used_live_time /(24.0*3600.0)
            print(f"New Amount of Events for {used_live_time_days} days ({used_live_time:.2f} s): {len(subset['time'])}")
            subset.meta['length'] = length

            # Save result
            write_table_hdf5(subset, output_file, path='timeseries', overwrite=True, serialize_meta = True)
            print(f"Written subset to {output_file} using {used_live_time_days:.4f} days of live time")

            # Create new GTI array corresponding to the extracted subset
            gti_output_file = os.path.join(
                data['output_dir'],
                f"{file_name}_{used_live_time_days:.2f}d_subset_gtis"
            )
            print(f"Used GTIs: {used_gtis}")
            saveGTIs(used_gtis, gti_output_file)
            print(f"Written GTI for subset to {gti_output_file}")





if __name__ == "__main__":
    main()