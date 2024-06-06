""" Performs barycentric and binary corrections on TimeSeries object generated from KM3NeT data.

Usage: CorrectEventListKM3NeT.py -i INPUT_FILES... [--wildcard] -o OUTPUT_DIR -s SOURCE_SPECS_FILE

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output directory
  -s --source SOURCE_SPECS_FILE          Hdf5 file containing information about the source of interest (ra, dec, P_orb, ...)     
"""
#python3 psr/src/scripts/CorrectEventListKM3NeT.py -i '/home/hpc/capn/capn107h/software/eventlistTestOutput/mcv*' -o correctedeventlistTestOutput/ -s hdf5SourceFiles/Vela_X-1.h5

from docopt import docopt
import os, glob
import h5py as h5py
import re
import numpy as np
from astropy.time import Time, TimeDelta
from astropy.timeseries import BinnedTimeSeries
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import EarthLocation, SkyCoord, Angle
from astropy.table import Table



# PLENS Imports
import plens.TimeSeries as TS
import plens.EventList as EL
from plens.TimeSeries import orbit_cor_deeter
import plens.antares_hdf5
import plens.antares_hdf5 as antares_hdf5




def main():
    # Getting Key-Argument-Pairs that are passed to the script
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

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])


    
    
    # Reading data files with source information
    with h5py.File(data['source'], 'r') as f:
        # Source name
        source_name = f['source_name'][()]
        print("Source Name:", source_name)

        # SkyCoord information
        ra = f['skycoord/ra'][()]
        dec = f['skycoord/dec'][()]
        skycoord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        print("SkyCoord:", skycoord)

        Porb = f['P_orb'][()]
        print("Orbital Period:", Porb, "days")
        
        Tpi2 = f['T_pi2'][()]
        axsini = f['axsini'][()]
        e = f['e'][()]
        omega = f['omega'][()]

    for file in input_files:
        # Fetching filename for usage in output filename 
        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]

        output_file = data['output_dir'] + file_name + '_corrected.hdf5'
        print(output_file, os.path.exists(output_file))

        # Read TimeSeries file and apply corrections using timing and source information
        with h5py.File(file) as input_file:
            # Read non corrected event list
            TimeSeries = EL.readEventList(input_file)
            print("Raw TimeSeries:")
            print(TimeSeries)
            # Barycentric correction
            TimeSeries = EL.barycentric_correction(TimeSeries, skycoord)
            print("After Baryocentric Correction:")
            print(TimeSeries)
            # Binary correction
            time_corr = orbit_cor_deeter(TimeSeries.time.value, 
                                            (float(Porb)*u.d).to_value(u.s), 
                                            float(axsini), 
                                            float(e), 
                                            Angle(float(omega), u.deg).radian - np.pi/2, 
                                            Time(float(Tpi2) + float(Porb)/2, format='jd').unix
                                        )
            # Creating new table with corrected eventlist
            TimeSeries = Table([time_corr], names=['time'])
            print("After Binary Corrections:")
            print(TimeSeries)

        # Store Corrected TimeSeries in HDF5-File
        with h5py.File(output_file, 'w') as output:   
            #EL.saveEventList(TimeSeries, output)
            TimeSeries.write(output, format='hdf5', overwrite=True, serialize_meta=True)
    
if __name__ == "__main__":
    main()
    