""" Epochfold corrected Eventlists and calculate the chi2 statistic for the testperiods.

Usage: FoldEventlist.py -i INPUT_FILES -o OUTPUT_DIR --filepattern=<filepattern> [--frequency=<frequency>] [--number_of_testf=<number_of_testf>] [--nbin=<nbin>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output directory
     --frequency=<float>                 Principle frequency around which an interval for the testfrequencies will                                          be chosen. [default: 10.]
     --number_of_testf=<float>             Number of testfrequencies to test around the principle frequency.                                                  [default: 200]
     --nbin=<int>                        Number of bins in the folded profile [default: 32]

"""
# python3 FoldEventlist.py -i './data/Antares_*_eventlist_*.hdf5' -o './data/' --filepattern 'Antares_(\d*)_eventlist_(\d*).hdf5'
import os, glob
from docopt import docopt
import re
import h5py
import numpy as np
import plens.TimeSeries as TS
from epochfolding.stingray_epochfolding import epochfolding_scan#, plot_efstat

def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    input_files = glob.glob(data['input_files'])
    input_files.sort()
    
    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])
        
    # Construct Eventlist for each run from corrected TimeSeries
    for file in input_files:
        #split = re.split('Antares_(\d*)_eventlist_(\d*).hdf5', file)
        split = re.split(data['filepattern'], file)
        run_number = split[1]
        split_number = split[2]
        print('Processing Run Nr.: ' + str(run_number) + ', split: ' + str(split_number), end='\n')
        
        output = data['output_dir'] + 'Antares_' + run_number + '_chi2_' + split_number + '.hdf5'
        
        if os.path.exists(output + '.hdf5'):
            continue
            
        with h5py.File(file) as input_file:
            
            _, effreq, efstat = epochfolding_scan(input_file, frequencies=float(data['frequency']), nbin=int(data['nbin']), oversampling=10, obs_length=100, number_testf=float(data['number_of_testf']), plot=False, save=True, format='hdf5', output=output)


if __name__ == "__main__":
    main()