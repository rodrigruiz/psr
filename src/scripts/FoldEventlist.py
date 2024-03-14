""" Epochfold corrected Eventlists and and save folded profiles.

Usage: FoldEventlist.py -i INPUT_FILES -o OUTPUT_DIR --filepattern=<filepattern> [--frequency=<frequency>] [--number_of_testf=<number_of_testf>] [--nbin=<nbin>] [--df=<float>] [--expocorr --gtis=<gtis>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output directory
     --frequency=<float>                 Principle frequency around which an interval for the testfrequencies will be chosen. [default: 10.]
     --number_of_testf=<int>             Number of testfrequencies to test around the principle frequency. [default: 200]
     --df=<float>                        Resolution of testfrequencies. [default: 1e-3]
     --nbin=<int>                        Number of bins in the folded profile. [default: 32]

"""
# python3 FoldEventlist.py -i './data/Antares_*_eventlist_*.hdf5' -o './data/' --filepattern 'Antares_(\d*)_eventlist_(\d*).hdf5'
import os, glob
from docopt import docopt
import re, fnmatch
import h5py
import numpy as np
import plens.TimeSeries as TS
import plens.EventList as EL
from epochfolding.stingray_epochfolding import epochfolding_single, get_testfrequencies#, plot_efstat
from epochfolding.gtis import loadGTIs

def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    input_files = glob.glob(data['input_files'])
    input_files.sort()
    #print(input_files)
    
    gti_files = glob.glob(data['gtis'])
    gti_files.sort()
    
    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])
        
    # Construct Eventlist for each run from corrected TimeSeries
    for file in input_files:
        #split = re.split('Antares_(\d*)_eventlist_(\d*).hdf5', file)
        split = re.split(data['filepattern'], file)
        run_number = split[1]
        split_number = split[2]
        print('Processing Run Nr.: ' + str(run_number) + ', split: ' + str(split_number), end='\n')
        
        gti_file = fnmatch.filter(gti_files, '*'+run_number+'*')
        print(gti_file)
        gtis = loadGTIs(gti_file[0])
        print(gtis)
        
        output = 'Antares_' + run_number + '_folded_' + split_number
        
        if os.path.exists(output + '.hdf5'):
            continue
            
        frequencies = get_testfrequencies(float(data['frequency']), int(data['number_of_testf']), float(data['df']))
        
        with h5py.File(file) as input_file:
            print(input_file.keys())
            epochfolding_single(input_file, 
                                frequencies, 
                                nbin=int(data['nbin']),
                                expocorr=data['expocorr'], 
                                gti=gtis, 
                                output=output, plot=False, save=True, 
                                format='hdf5', outputdir=data['output_dir']
                               )
            
            
if __name__ == "__main__":
    main()