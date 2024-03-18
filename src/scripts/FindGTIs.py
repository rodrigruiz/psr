""" Check for gaps in the data. Determine good time intervals.

Usage: FindGTIs.py <INPUT_FILES>... [--wildcard] -o OUTPUT_DIR --filepattern=<filepattern> [--count --outputfile=<outputfile>] [--combine]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output directory
     --find                              Whether to find gtis in Antares runs.
     --combine                           Whether to save gtis in one file.
  
"""
# 

from docopt import docopt
import os, glob
import h5py
import re
import pickle
import numpy as np
import plens.TimeSeries as TS
from astropy.table import Table
from epochfolding.gtis import findGTIs, saveGTIs

def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    if data['wildcard']:
        input_files = glob.glob(data['<INPUT_FILES>'][0])
    else:
        input_files = data['<INPUT_FILES>']
    input_files.sort()
    
    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])
    
    if data['count']:
        Ngtis = Table(names=('run_number', 'N_bad_data_rows'), dtype=(str, int))
    
    if data['combine']:
        gtis = []
        
    for file in input_files:
        split = re.split(data['filepattern'], file)
        run_number = split[1]
        
        #print(split)
        if split[2] != '':
            split_number = split[2]
            print('Processing Run Nr.: ' + str(run_number) + ', split: ' + str(split_number), end='\n')
            
            if data['combine']:
                output = data['output_dir'] + 'Antares_gtis'
            else:
                output = data['output_dir'] + 'Antares_' + run_number + '_gtis_' + split_number
            
        else:
            print('Processing Run Nr.: ' + str(run_number), end='\n')
            
            if data['combine']:
                output = data['output_dir'] + 'Antares_gtis'
            else:
                output = data['output_dir'] + 'Antares_' + run_number + '_gtis'
        
            
        with h5py.File(file) as input_file:
            
            ts, timeslice_duration = TS.readTimeSeries(input_file)
            
            if data['combine']:
                gtis += findGTIs(ts['rateOn'].value, ts['time_bin_start'].value)
            else:
                gtis = findGTIs(ts['rateOn'].value, ts['time_bin_start'].value)
                saveGTIs(gtis, output)
            
            if data['count']:
                Ngtis.add_row((run_number, len(ts[ ts['rateOn'].value == 0. ])))

    if data['combine']:
        gtis.sort()
        saveGTIs(gtis, output)
                                 
    if data['count']:
        Ngtis.write(data['output_dir'] + data['outputfile'], format='ascii', overwrite=True)
    
if __name__ == "__main__":
    main()
    