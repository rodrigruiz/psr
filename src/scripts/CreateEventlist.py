""" Create Eventlist from TimeSeries

Usage: CreateEventlist.py -i INPUT_FILES -o OUTPUT_DIR --filepattern=<filepattern> --scaling=<scaling> --overwrite=<BOOL>

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output directory
     --scaling=<scaling>                 Whether to scale the ANTARES rates. 
                                         The actual rate will be divided by this number.

"""
# python3 CreateEventlist.py -i './data/Antares_*_total_rates_*_corrected.hdf5' -o './data/' --scaling 1e10 --filepattern 'Antares_(\d*)_total_rates_(\d*)_corrected.hdf5'

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
    #print(input_files)
    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])
    
    # Construct Eventlist for each run from corrected TimeSeries
    for file in input_files:
        #split = re.split('Antares_(\d*)_total_rates_(\d*)_corrected.hdf5', file)
        #split = re.split('Antares_(\d*)_total_rates_(\d*).hdf5', file)
        split = re.split(data['filepattern'], file)
        run_number = split[1]
        split_number = split[2]
        print('Processing Run Nr.: ' + str(run_number) + ', split: ' + str(split_number), end='\n')
        
        output = data['output_dir'] + 'Antares_' + run_number + '_eventlist_' + split_number
        #print(output, os.path.exists(output + '.hdf5'))
        
        #if data['overwrite'] is False:
        if os.path.exists(output + '.hdf5'):
            continue
            
        with h5py.File(file) as input_file:
            ts, timeslice_duration = TS.readTimeSeries(input_file)  
            create_eventlist(ts=ts, output=output, 
                         time_column_name='time_bin_start', 
                         data_column_name='rateOn', 
                         scaling=float(data['scaling']),
                         format='hdf5')
            

if __name__ == "__main__":
    main()
    