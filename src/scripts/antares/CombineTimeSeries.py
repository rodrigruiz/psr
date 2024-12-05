""" Combines Time Series

Usage: CombineTimeSeries.py -i INPUT_FILES -o OUTPUT_FILE

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_file OUTPUT_FILE           Output file

"""
# python3 CombineTimeSeries.py -i "/home/saturn/capn/mppi19/data/ANTARES/rates/hdf5/Antares_0522*_total_rates.hdf5" -o 'CombinedTimeSeries_Test.hdf5'
# python3 CombineTimeSeries.py -i "/home/saturn/capn/mppi19/data/ANTARES/rates/hdf5/Antares_*_total_rates.hdf5" -o 'CombinedTimeSeries_MEff_MA.hdf5'
from docopt import docopt
import glob
import re
import h5py as h5py

from astropy.table import vstack

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

    rates = list()
    
    # Construct TimeSeries for each run
    for file in input_files:
        with h5py.File(file) as h5_file:
            run_number = re.split('Antares_(\d*)_total_rates.hdf5', file)[1]
            print('Processing Run Nr.: ' + str(run_number), end='\r')
            rates.append(plens.TimeSeries.get_TimeSeries(h5_file))
            timeslice_duration = plens.antares_hdf5.get_timeslice_duration(h5_file)
        
    # Cobine TimeSeries
    TimeSeries = vstack(rates)

    # Write Output File
    with h5py.File(data['output_file'], 'w') as output_file:   
    
        plens.TimeSeries.saveTimeSeries(TimeSeries, timeslice_duration, output_file)

 
            
if __name__ == "__main__":
    main()
