""" Performs corrections on TimeSeries.

Usage: CorrectTimeSeries.py -i INPUT_FILES -o OUTPUT_DIR --correction=<correction> [--rajd RIGHT_ASCENSION --decjd DECLINATION]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output file
     --correction=<correction>           Which correction to perform on the TimeSeries. Options: 'bary', 'fill', 'all' [default: all]
                                         If 'bary' or 'all' is selected RA and Dec should also be specified.
     --rajd RIGHT_ASCENSION              Right ascension (J2000) (degrees)
     --decjd DECLINATION                 Declination (J2000) (degrees)
     
"""
#python3 CorrectTimeSeries.py -i '../stingray/antares/data/Antares_051870_total_rates_*[!combined].hdf5' -o '../stingray/antares/data/' --correction bary --rajd 02h43m40.4252869512s --decjd +61d26m03.757456824s
# python3 CorrectTimeSeries.py -i 'CombinedTimeSeries_Test.hdf5' -o 'CorrectedTimeSeries_Test.hdf5' --rajd '256.06166667' --decjd '-60.28166667'
# python3 CorrectTimeSeries.py -i "CombinedTimeSeries_MEff_MA.hdf5" -o "TimeSeries_J1704-6016.hdf5" --rajd "256.06166667" --decjd "-60.28166667"

from docopt import docopt
import os, glob
import h5py as h5py
import re

from astropy.time import Time, TimeDelta
from astropy.timeseries import BinnedTimeSeries
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import EarthLocation


# PLENS Imports
import plens.TimeSeries as TS
import plens.antares_hdf5
import plens.antares_hdf5 as antares_hdf5




def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    if data['rajd'] != None and data['decjd'] != None:
        skycoord = coord.SkyCoord(ra=data['rajd'], dec=data['decjd'], unit=(u.deg, u.deg))
    else:
        raise ValueError(
            "please specify RA and Dec of the object!")
    
    input_files = glob.glob(data['input_files'])
    input_files.sort()
    
    # Construct TimeSeries for each run
    for file in input_files:
        # Read TimeSeries
        with h5py.File(file) as input_file:
            split = re.split('Antares_(\d*)_total_rates_(\d*).hdf5', file)
            run_number = split[1]
            split_number = split[2]
            
            print('Processing Run Nr.: ' + str(run_number) + ', split: ' + str(split_number), end='\r')
            output_file = data['output_dir'] + 'Antares_' + run_number + '_total_rates_' + split_number + '_corrected.hdf5'
            
            TimeSeries, timeslice_duration = TS.readTimeSeries(input_file)  
            
            
    
            # Perform Barycentric Correction
            if data['correction'] == 'bary':
           
                #print('Perform Barycentric Correction')
                TimeSeries = TS.barycentric_correction(TimeSeries, timeslice_duration, skycoord)
                #print('Barycentric Correction Successful')
            
            
            elif data['correction'] == 'fill':
                # Fill Gaps in TimeSeries
                TimeSeries = TS.fill_TimeSeries_gaps(TimeSeries, timeslice_duration)
        
            elif data['correction'] == 'all':
        
                #print('Perform Barycentric Correction')
                TimeSeries = TS.barycentric_correction(TimeSeries, timeslice_duration, skycoord)
                #print('Barycentric Correction Successful')
        
                # Fill Gaps in TimeSeries
                TimeSeries = TS.fill_TimeSeries_gaps(TimeSeries, timeslice_duration)
    
            else:
                raise ValueError(
                    "correction was not specified or the correction was not recognised.")
    
        # Store Corrected TimeSeries again in HDF5-File
        with h5py.File(output_file, 'w') as output:  
            
            TS.saveTimeSeries(TimeSeries, timeslice_duration, output)
        
    
if __name__ == "__main__":
    main()
    