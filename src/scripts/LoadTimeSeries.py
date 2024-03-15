""" Load Time Series

Usage: LoadTimeSeries.py -i INPUT_FILES -o OUTPUT_DIR [--blind]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output directory
     --blind                             Whether to blind the data.

"""
#/home/wecapstor3/capn/mppi19/ANTARES/data/rates/hdf5/
#python3 antares_loadtimeseries.py -i '/home/wecapstor3/capn/mppi19/ANTARES/data/rates/hdf5/Antares_05*0_total_rates.hdf5' -o './data/'

from docopt import docopt
import os, glob
import re
import h5py as h5py

# PLENS Imports
import plens.TimeSeries
import plens.antares_hdf5



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
        with h5py.File(file) as h5_file:
            run_number = re.split('Antares_(\d*)_total_rates.hdf5', file)[1]
            print('Processing Run Nr.: ' + str(run_number), end='\n')
            
            output_file = data['output_dir'] + 'Antares_' + run_number + '_total_rates_combined.hdf5'
            
            TimeSeries = plens.TimeSeries.get_TimeSeries(h5_file)
            timeslice_duration = plens.antares_hdf5.get_timeslice_duration(h5_file)
            
            # blind data
            if data['blind']:
                plens.TimeSeries.shuffleTimeSeries(TimeSeries, 'rateOn')
                
            # Write Output File
            with h5py.File(output_file, 'w') as output_file:   
                plens.TimeSeries.saveTimeSeries(TimeSeries, timeslice_duration, output_file)

 
            
if __name__ == "__main__":
    main()
