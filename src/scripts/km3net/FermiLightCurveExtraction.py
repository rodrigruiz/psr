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

    # Initialize combined_time and combined_energy as empty numpy arrays
    combined_time = np.array([])
    combined_energy = np.array([])
    combined_ra = np.array([])
    combined_dec = np.array([])
    combined_gti_start = np.array([])
    combined_gti_stop = np.array([])

    for file in input_files:
        with fits.open(file) as hdul:
                data = hdul[1].data  # Assuming light curve data is in the first extension
                ext1 = hdul[2].data
                #print(data.__dict__)
                #print(ext1.__dict__)
                time = data['TIME']
                energy = data['ENERGY']
                ra = data['RA']
                dec = data['DEC']

                gti_start = ext1['START']
                gti_stop = ext1['STOP']
                print(f"File: {file} - Time entries: {len(time)}")

                # Concatenate the data from each file
                combined_time = np.concatenate((combined_time, time))
                combined_energy = np.concatenate((combined_energy, energy))
                combined_ra = np.concatenate((combined_ra, ra))
                combined_dec = np.concatenate((combined_dec, dec))
                combined_gti_start = np.concatenate((combined_gti_start, gti_start))
                combined_gti_stop = np.concatenate((combined_gti_stop, gti_stop))
                print(f"Combined time length: {len(combined_time)}")

    combined_unix_time = convert_mission_time_to_unix(combined_time)
    # Extract common prefix from input filenames
    common_prefix = os.path.commonprefix(input_files)
    # Remove any trailing non-alphanumeric characters from common prefix
    common_prefix = os.path.basename(common_prefix).rstrip("_-.")
    if not common_prefix:
        common_prefix = "Test"
    output_file = os.path.join(output_dir, f"{common_prefix}_fermi_gammaray_lightcurve")


    # Create the output table
    event_list = Table([combined_unix_time, combined_energy, combined_ra, combined_dec],names=('time','energy','ra','dec'))
    gti_list = Table([combined_gti_start,combined_gti_stop],names=('gti_start','gti_stop'))


    print(event_list)
    print(gti_list)
    # Save the output table
    
    event_list.write(output_file + '.hdf5', format='hdf5', overwrite=True, serialize_meta=True)
    gti_list.write(output_file + '_gtilist.hdf5', format='hdf5', overwrite=True, serialize_meta=True)

    
if __name__ == "__main__":
    main()