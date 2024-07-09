""" Epochfold corrected KM3NeT Eventlists and and save folded profiles .

Usage: FoldEventListKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR [--frequency=<frequency>] [--number_of_testf=<number_of_testf>] [--nbin=<nbin>] [--df=<float>] [--expocorr --gtis=<gtis>] [--signal_strength=<signal_strength>]

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
#import stingray
def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]

    input_files = []
    for pattern in data['input_files']:
        input_files.extend(glob.glob(pattern))
    input_files.sort()

    if not input_files:
        print(f"No files matching pattern: {input_files}")
        return

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])
        
    # Construct Eventlist for each run from corrected TimeSeries
    for file in input_files:
        # Fetching filename for usage in output filename 
        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]
        #print(file_name)

        output_file = file_name + '_epochfolding'  
            
        frequencies = get_testfrequencies(float(data['frequency']), int(data['number_of_testf']), float(data['df']))
        
        with h5py.File(file) as input_file:
            #print(input_file.keys())
            epochfolding_single(input_file, 
                                frequencies, 
                                nbin=int(data['nbin']),
                                expocorr=data['expocorr'], 
                                #gti=gtis, 
                                output=output_file, plot=False, save=True, 
                                format='hdf5', outputdir=data['output_dir']
                               )
            
            
if __name__ == "__main__":
    main()
