""" Create Eventlist from TimeSeries

Usage: CreateEventlist.py -i INPUT_FILES -o OUTPUT_DIR --scaling=<scaling>

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output directory
     --scaling=<scaling>                 Whether to scale the ANTARES rates. 
                                         The actual rate will be divided by this number.

"""
# python3 CreateEventlist.py -i './data/Antares_*_total_rates_*_corrected.hdf5' -o './data/' --scaling 1e10
from docopt import docopt
import os, glob
import h5py
import re
import numpy as np
import plens.TimeSeries as TS
from epochfolding.generate_eventsfromtimeseries import create_eventlist

def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    input_files = glob.glob(data['input_files'])
    input_files.sort()
    
    # Construct Eventlist for each run from corrected TimeSeries
    for file in input_files:
        
        with h5py.File(file) as input_file:
            split = re.split('Antares_(\d*)_total_rates_(\d*)_corrected.hdf5', file)
            run_number = split[1]
            split_number = split[2]
        
            output = data['output_dir'] + 'Antares_' + run_number + '_eventlist_' + split_number
            
            ts, timeslice_duration = TS.readTimeSeries(input_file)
            
            create_eventlist(ts=ts, output=output, 
                         time_column_name='time_bin_start', 
                         data_column_name='rateOn', 
                         scaling=float(data['scaling']),
                         format='hdf5')

if __name__ == "__main__":
    main()
    