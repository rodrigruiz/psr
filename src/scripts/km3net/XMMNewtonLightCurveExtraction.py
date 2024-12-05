""" Loads and combines downloaded photon files from the Fermi LAT instrument and generates EventList for usage in the epoch folding analysis.

Usage: FermiLightCurveExtraction.py -i INPUT_FILES... -o OUTPUT_DIR

Options:
  -h --help                              Show this help message
  -i --input_files INPUT_FILES...        Input files
  -o --output_dir OUTPUT_DIR             Output directory  
"""

# python3 psr/src/scripts/CombineEventListsKM3NeT.py -i "/home/hpc/capn/capn107h/software/correctedeventlistTestOutput/*" -o combinedeventlistTestOutput/


from docopt import docopt
import os
import glob
import h5py
import numpy as np

from astropy.table import vstack, Table
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
from astropy.time import Time

from stingray.events import EventList
from stingray import Lightcurve

import os
import warnings

def convert_mission_time_to_unix(seconds_since_mission_start):
    """
    Convert Fermi MET (seconds since January 1, 2001) to Unix timestamp format.
    """
    # Define Fermi mission start time in MJD
    mjdref = 51910.00074287037  # Fermi mission reference in MJD
    # Convert Fermi mission reference to Unix time
    mission_start_unix = Time(mjdref, format='mjd', scale='tt').unix
    # Add seconds since mission start to mission_start_unix to get the Unix timestamps
    unix_timestamps = mission_start_unix + seconds_since_mission_start
    return unix_timestamps

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
    
    output_dir = data['output_dir']
    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    combined_time = np.array([])
    combined_rate = np.array([])
    combined_error = np.array([])

    for file in input_files:
        with fits.open(file) as hdul:
            data = hdul[1].data  # Assuming light curve data is in the first extension
            print(data.__dict__)
            time = data['TIME']
            rate = data['RATE']
            error = data['ERROR']
            print(f"File: {file} - Time entries: {len(time)} - Min Time: {np.min(time)} - Max Time: {np.max(time)}")
            
            # Concatenate the data from each file
            combined_time = np.concatenate((combined_time, time))
            combined_rate = np.concatenate((combined_rate, rate))
            combined_error = np.concatenate((combined_error, error))
            print(f"Combined time length: {len(combined_time)}")
            print(f"rate: {rate}")

    lc = Lightcurve(combined_time,combined_rate,combined_error, input_counts=False)
    events = EventList()
    events.simulate_times(lc)
    print(len(events.time))
    print(events.time)
    events_time = events.time
    combined_unix_time = convert_mission_time_to_unix(events_time)
    # Extract common prefix from input filenames
    common_prefix = os.path.commonprefix(input_files)
    # Remove any trailing non-alphanumeric characters from common prefix
    common_prefix = os.path.basename(common_prefix).rstrip("_-.")
    if not common_prefix:
        common_prefix = "Test"
    output_file = os.path.join(output_dir, f"{common_prefix}_xmmnewton_xray_lightcurve")

    print(combined_unix_time)


    # Create the output table
    event_list = Table([combined_unix_time],names='time')
    #gti_list = Table([combined_gti_start,combined_gti_stop],names=('gti_start','gti_stop'))


    print(event_list)
    #print(gti_list)
    # Save the output table
    
    event_list.write(output_file + '.hdf5', format='hdf5', overwrite=True, serialize_meta=True)
    #gti_list.write(output_file + '_gtilist.hdf5', format='hdf5', overwrite=True, serialize_meta=True)

    
if __name__ == "__main__":
    main()