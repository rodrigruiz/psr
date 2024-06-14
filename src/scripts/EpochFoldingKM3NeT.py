"""Apply stingray Epoch Folding Search to KM3NeT data.
Usage: EpochFoldingKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR [--frequency=<frequency>] [--number_of_testf=<number_of_testf>] [--nbin=<nbin>] [--df=<float>] 

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output file
     --frequency=<float>                 Principle frequency around which an interval for the testfrequencies will be chosen. [default: 10.]
     --number_of_testf=<int>             Number of testfrequencies to test around the principle frequency. [default: 200]
     --df=<float>                        Resolution of testfrequencies. [default: 1e-5]
     --nbin=<int>                        Number of bins in the folded profile. [default: 32]
"""

from docopt import docopt
import os, glob
import h5py as h5py
import re
import numpy as np
from matplotlib import pyplot as plt

# PLENS Imports
import plens.EventList as EL

from stingray.pulse.search import epoch_folding_search, z_n_search
from epochfolding.stingray_epochfolding import savehdf5, get_testfrequencies

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

    for file in input_files:    
        # Fetching filename for usage in output filename 
        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]
        #print(file_name)

        output_plot = data['output_dir'] + file_name + '_epochfolding_resultplot.png'
        output_file = data['output_dir'] + file_name + '_epochfolding_results.hdf5'

        with h5py.File(file) as input_file:
                
            EventList = EL.readEventList(input_file)
            print(EventList)
            frequencies = get_testfrequencies(float(data['frequency']), int(data['number_of_testf']), float(data['df']))
            print(frequencies)
            freq, efstat = epoch_folding_search(np.array(EventList['time'].value), frequencies, nbin=int(data['nbin']))

            # ---- PLOTTING --------
            plt.figure()
            plt.plot(freq, efstat, label='EF statistics')
            plt.axhline(int(data['nbin']) - 1, ls='--', lw=3, color='k', label='n - 1')
            plt.axvline(float(data['frequency']), lw=3, alpha=0.5, color='r', label='Correct frequency')
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('EF Statistics')
            _ = plt.legend()
            plt.savefig(output_plot)

            with h5py.File(output_file, 'w') as out:
                savehdf5(freq, efstat, out)

if __name__ == "__main__":
    main()
