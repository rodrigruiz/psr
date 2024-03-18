""" Split Time Series

Usage: SplitTimeSeries.py -i INPUT_FILES -o OUTPUT_DIR --bins_per_file=<bins_per_files>

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output directory
     --bins_per_file=<bins_per_files>    Bins per file
                                         Splits the original timeseries in total_bins/bins_per_file files.

"""
# python3 SplitTimeSeries.py -i './data/Antares_05*_total_rates_combined.hdf5' -o './data/' --bins_per_file 100

import h5py
import numpy as np
from docopt import docopt
import os, glob
import re
import plens.TimeSeries as TS
from epochfolding.generate_split_timeseries import split_timeseries

def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    
    input_files = glob.glob(data['input_files'])
    input_files.sort()

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])
        
    # Construct TimeSeries for each run
    for file in input_files:
        
        run_number = re.split('Antares_(\d*)_total_rates_combined.hdf5', file)[1]
        print('Processing Run Nr.: ' + str(run_number), end='\n')
        
        output_file = data['output_dir'] + 'Antares_' + run_number + '_total_rates'
            
        if os.path.exists(output_file + '.hdf5'):
            continue
                
        with h5py.File(file) as h5_file:
    
            ts, timeslice_duration = TS.readTimeSeries(h5_file)
            
            split_timeseries(ts, len(ts), np.ceil(len(ts)/int(data['bins_per_file'])), output_file, format='hdf5', antares=True, timeslice_duration=timeslice_duration)

if __name__ == "__main__":
    main()
    