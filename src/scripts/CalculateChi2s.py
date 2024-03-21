""" Epochfold corrected Eventlists and calculate the chi2 statistic for the testperiods.

Usage: CalculateChi2s.py <INPUT_FILES>... [--wildcard] -o OUTPUT_DIR --filepattern=<filepattern> [--frequency=<frequency>] [--number_of_testf=<number_of_testf>] [--nbin=<nbin>] [--expocorr --gti_files=<gti_files>]

Options:
  -h --help                              Help
  -i --input_files <INPUT_FILES>           Input files
  -o --output_dir OUTPUT_DIR             Output directory
     --frequency=<float>                 Principle frequency around which an interval for the testfrequencies will be chosen. [default: 10.]
     --number_of_testf=<float>             Number of testfrequencies to test around the principle frequency. [default: 200]
     --nbin=<int>                        Number of bins in the folded profile [default: 32]
     --expocorr                          Whether to perform exposire correction. [default: False]
     --gti_files=<GTI_FILES>                      GTI files. works with only one gti file or a wildcard expression.

"""
# python3 FoldEventlist.py -i './data/Antares_*_eventlist_*.hdf5' -o './data/' --filepattern 'Antares_(\d*)_eventlist_(\d*).hdf5'
import os, glob
from docopt import docopt
import re, fnmatch
import h5py
import numpy as np
import plens.TimeSeries as TS
import plens.EventList as EL
from epochfolding.stingray_epochfolding import epochfolding_scan
from epochfolding.gtis import loadGTIs

def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    print(data)
    if data['wildcard']:
        #input_files = glob.glob(data['input_files'][0])
        input_files = glob.glob(data['<INPUT_FILES>'][0])
        #gti_files = glob.glob(data['gti_files'][0])
    else:
        input_files = data['<INPUT_FILES>']
        #gti_files = data['gti_files']
    input_files.sort()
    
    if data['expocorr']:
        gti_files = glob.glob(data['gti_files'])
        print(gti_files)
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
        
        #print(gti_file)
        if data['expocorr']:
            if len(gti_files) > 1:
                gti_file = fnmatch.filter(gti_files, '*'+run_number+'*')
            #print(gti_file)
        #else:
                gtis = loadGTIs(gti_file[0])
            #print(gtis)
            else:
                gtis = loadGTIs(gti_files[0])
        else:
            gtis = None
            
        output = data['output_dir'] + 'Antares_' + run_number + '_chi2_' + split_number + '.hdf5'
        
        if os.path.exists(output + '.hdf5'):
            continue
            
        with h5py.File(file) as input_file:
            #print(input_file.keys(), input_file['timeseries'].keys())
            #print(input_file.keys(), input_file['__astropy_table__'])
            timeseries = EL.readEventList(input_file)
            
            #print(re.search('*'+run_number+'*', gti_files))
            #search = [re.search(str(run_number), gti_file) for gti_file in gti_files]
            #print(search)
            
           
            # calculating statistics for chi2 histogram
            _, effreq, efstat = epochfolding_scan(input_file, frequencies=float(data['frequency']), nbin=int(data['nbin']), oversampling=10, obs_length=100, number_testf=float(data['number_of_testf']), expocorr=data['expocorr'], gti=gtis, plot=False, save=True, format='hdf5', output=output)


if __name__ == "__main__":
    main()