""" Epochfold corrected Eventlists and and save folded profiles.

Usage: EFAddedProfiles.py -i INPUT_DIR [--wildcard] -o OUTPUT_DIR --filepattern=<filepattern> [--frequency=<frequency>] [--number_of_testf=<number_of_testf>] [--nbin=<nbin>] [--df=<float>] [--signal_strength=<signal_strength>]

Options:
  -h --help                              Help
  -i --input_dir INPUT_DIR           Input directory
  -o --output_dir OUTPUT_DIR             Output directory
     --frequency=<float>                 Principle frequency around which an interval for the testfrequencies will be chosen. [default: 10.]
     --number_of_testf=<int>             Number of testfrequencies to test around the principle frequency. [default: 200]
     --df=<float>                        Resolution of testfrequencies. [default: 1e-3]
     --nbin=<int>                        Number of bins in the folded profile. [default: 32]

"""
# python3 FoldEventlist.py -i './data/Antares_*_eventlist_*.hdf5' -o './data/' --filepattern 'Antares_(\d*)_eventlist_(\d*).hdf5'
import os, glob
from docopt import docopt
import re
import h5py
import numpy as np
import plens.TimeSeries as TS
import plens.EventList as EL
#from stingray.events import EventList
from epochfolding.ef_addedprofiles import epoch_folding_search
from epochfolding.stingray_epochfolding import savehdf5, get_testfrequencies

def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    #gti_files = glob.glob(data['gtis'])
    #gti_files.sort()
    
    """
    if data['wildcard']:
        #input_files = glob.glob(data['input_files'][0])
        input_files = glob.glob(data['<INPUT_FILES>'][0])
        #gti_files = glob.glob(data['gti_files'][0])
    else:
        input_files = data['<INPUT_FILES>']
        #gti_files = data['gti_files']
    input_files.sort()
    """
    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])
    
    #gti_file = fnmatch.filter(gti_files, '*'+run_number+'*')
    #print(gti_file)
    #gtis = loadGTIs(gti_file[0])
    #print(gtis)
    #for file in input_files:
    #    split = re.split(data['filepattern'], file)
    #    run_number = split[1]
    #    split_number = split[2]
    #    print('Processing Run Nr.: ' + str(run_number) + ', split: ' + str(split_number), end='\n')
    output = data['output_dir'] + 'Antares_chi2_added_' + data['signal_strength'] + '.hdf5'
        
    #if os.path.exists(output + '.hdf5')
        
    frequencies = get_testfrequencies(float(data['frequency']), int(data['number_of_testf']), float(data['df']))
    #frequencies = [9.9, 9.901, 9.902]    
        #for p in frequencies:
    #with h5py.File() as input_file:
    #    events = EventList().read(input_file, 'hdf5')
    
    # dummy variable
    times = np.arange(1.28408629e+09, 1.28408630e+09, 0.1)   
    
    # all splits for a certain period
    #files_to_add = glob.glob()
    effreq, efstat = epoch_folding_search(times, 
                                 frequencies, 
                                 data['input_dir'], 
                                 data['filepattern'],
                                 data['output_dir'],
                                 nbin=int(data['nbin']), 
                                 segment_size=np.inf, 
                                 #expocorr=data['expocorr'], 
                                 #gti=gtis,
                                 expocorr=False,
                                 gti=None,
                                 weights=1, 
                                 fdots=0
                                )
    
    with h5py.File(output, 'w') as out:
        savehdf5(effreq, efstat, out)
            
if __name__ == "__main__":
    main()