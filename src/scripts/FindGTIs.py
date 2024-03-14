""" Check for gaps in the data. Determine good time intervals.

Usage: FindGTIs.py -i INPUT_FILES -o OUTPUT_DIR --filepattern=<filepattern> [--count --outputfile=<outputfile>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output directory
     --find                              Whether to find gtis in Antares runs.
  
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

def findGTIs(a, time, outputfile):
    """Finds 0. entries in data array and returns time intervals of good data.
    
    Parameters
    ----------
        a : np.array
            Data to test.
        
        time : np.array
            Times for the data.
            
    Returns
    -------
        gtis : np.array
            Good Time Intervals.
        
    """
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    rate_zero = np.where(absdiff == 1)[0].reshape(-1, 2)
    #print(rate_zero)
    # correct 0 column entries for gtis
    for index, i in enumerate(rate_zero[:,0]):
        rate_zero[index][0] = i - 1
    #print(rate_zero)
    
    # 0th and last element
    indices = np.concatenate(([0], rate_zero.flatten(), [len(a)-1])).reshape(-1, 2)
    
    # check if first bin is bti
    if (indices[0] < 0.).any():
        indices = np.delete(indices, 0, 0)
    if (indices[-1] >= len(a)).any():
        indices = np.delete(indices, -1, 0)

    return list((map(tuple, time[indices])))

def saveGTIs(gtis, outputfile):
    """Save GTIs to pickle file.
    
    Parameters
    ----------
        gtis : np.array
            Good Time Intervals.
        
        outputfile : str
            Filename/Filpath to save GTIs to. Without extension.
            
    """
    
    with open(outputfile + '.pkl', 'wb') as f:
        pickle.dump(gtis, f)
        
    #with open('testpickle.pkl', 'rb') as f:
    #gtis_file = pickle.load(f)
    
    return

def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    input_files = glob.glob(data['input_files'])
    input_files.sort()
    
    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])
    
    if data['count']:
        Ngtis = Table(names=('run_number', 'N_bad_data_rows'), dtype=(str, int))

    for file in input_files:
        split = re.split(data['filepattern'], file)
        run_number = split[1]
        
        #print(split)
        if split[2] != '':
            split_number = split[2]
            print('Processing Run Nr.: ' + str(run_number) + ', split: ' + str(split_number), end='\n')
            
            output = data['output_dir'] + 'Antares_' + run_number + '_gtis_' + split_number + '.pkl'
            
        else:
            print('Processing Run Nr.: ' + str(run_number), end='\n')
            
            output = data['output_dir'] + 'Antares_' + run_number + '_gtis.pkl'
        
            
        with h5py.File(file) as input_file:
            
            ts, timeslice_duration = TS.readTimeSeries(input_file)
            
            gtis = findGTIs(ts['rateOn'].value, ts['time_bin_start'].value, output)
            
            saveGTIs(gtis, output)
            
            if data['count']:
                Ngtis.add_row((run_number, len(ts[ ts['rateOn'].value == 0. ])))

    
    if data['count']:
        Ngtis.write(data['output_dir'] + data['outputfile'], format='ascii', overwrite=True)
    
if __name__ == "__main__":
    main()
    