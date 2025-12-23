"""Apply stingray Epoch Folding Search to KM3NeT data.
Usage: EpochFoldingKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR [--gti_files=<gti_files>...] [--frequency=<frequency>] [--number_of_testf=<number_of_testf>] [--nbin=<nbin>] [--df=<float>] [--ratio=<ratio>] [--iteration=<iteration>] [--segment_size=<segment_size>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output file
     --gti_files=<gti_files>...          Optional GTI files
     --frequency=<float>                 Principle frequency around which an interval for the testfrequencies will be chosen. [default: 10.]
     --number_of_testf=<int>             Number of testfrequencies to test around the principle frequency. [default: 200]
     --df=<float>                        Resolution of testfrequencies. [default: 1e-5]
     --nbin=<int>                        Number of bins in the folded profile. [default: 32]
     --ratio=<float>                     Signal to Noise ratio [default: 0.3]
     --iteration=<int>                   Nr of current iteration [default: 0]
     --segment_size=<float>              Length of the segments to be averaged in the periodogram [default: 5000]
"""

from docopt import docopt
import os, glob
import h5py as h5py
import re
import numpy as np
from matplotlib import pyplot as plt

import stingray
import inspect

from astropy.io.misc.hdf5 import read_table_hdf5
from astropy.timeseries import TimeSeries


from stingray.pulse.search import epoch_folding_search, z_n_search, search_best_peaks
from epochfolding.stingray_epochfolding import savehdf5, get_testfrequencies
from epochfolding.gtis import loadGTIs

def readEventList(file):
    """Read EventList from HDF5-file constucted by CreateEventlist. 
    
    Parameters
    ----------
    file : h5py.File
        Input file to read EventList.
    
    Returns
    -------
    timeseries : astropy.timeseries.BinnedTimeSeries
        The BinnedTimeSeries has four columns: 'time_bin_start'
    
    """

    timeseries = TimeSeries.read(file, format='hdf5', time_column='time', time_format='unix')

    
    return timeseries

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
    
    gti_files = data.get('gti_files', [])
    if gti_files is None:
        gti_files = []
    else:
        gti_files = [file for file in gti_files]

    gti_files.sort()

    # Validate GTI file count
    if len(gti_files) == 1 and len(input_files) > 1:
        print("Warning: Only one GTI file provided for multiple input files. It will be applied to all input files.")
    elif len(gti_files) not in [0, 1, len(input_files)]:
        print("Error: Number of GTI files must be either 0, 1, or equal to the number of input files.")
        return

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    expocorr = False

    for idx, file in enumerate(input_files):
        # Determine GTI for the current input file
        current_gti = None
        if gti_files:
            expocorr = True
            gti_file = gti_files[0] if len(gti_files) == 1 else gti_files[idx]

            #gti_table = read_table_hdf5(gti_file)
            print(gti_file)
            gti_table = loadGTIs(gti_file)
            #gti_table = np.array(gti_table)
            gti_start = gti_table[:,0]
            gti_stop = gti_table[:,1]
            current_gti = np.array([gti_start, gti_stop]).T
            print(f"GTIs: {current_gti}")
            

        # Fetching filename for usage in output filename
        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]

        #output_plot = os.path.join(data['output_dir'], f"{file_name}_SNR_{data['ratio']}_I{data['iteration'].zfill(4)}_TestFrequency_{data['frequency']}_epochfolding_resultplot.png")
        #output_file = os.path.join(data['output_dir'], f"{file_name}_SNR_{data['ratio']}_I{data['iteration'].zfill(4)}_TestFrequency_{data['frequency']}_epochfolding_results.hdf5")
        output_file = os.path.join(data['output_dir'], f"{file_name}_r{data['ratio']}_I{data['iteration'].zfill(4)}_{data['frequency']}Hz_ef.hdf5")
        output_plot = os.path.join(data['output_dir'], f"{file_name}_r{data['ratio']}_I{data['iteration'].zfill(4)}_{data['frequency']}Hz_efplot.png")
        #output_plot = os.path.join(data['output_dir'], f"SNR_{data['ratio']}_I{data['iteration'].zfill(4)}_TestFrequency_{data['frequency']}_epochfolding_resultplot.png")
        #output_file = os.path.join(data['output_dir'], f"SNR_{data['ratio']}_I{data['iteration'].zfill(4)}_TestFrequency_{data['frequency']}_epochfolding_results.hdf5")

        with h5py.File(file, 'r') as input_file:
            EventList = readEventList(input_file)
            print(f"Times: {np.array(EventList['time'].value)}")
            frequencies = get_testfrequencies(float(data['frequency']), int(data['number_of_testf']), float(data['df']))
            #print(f"First 10 Current GTIs: {current_gti[:10]}")
            #print(f"Last 10 Current GTIs: {current_gti[-10:]}")
            print(f"First 10 Times: {[f'{x:.3f}' for x in np.array(EventList['time'].value)[:10]]}")
            print(f"Last 10 Times: {[f'{x:.3f}' for x in np.array(EventList['time'].value)[-10:]]}")
            print(f"Number of available events: {len(np.array(EventList['time'].value))}")
            print(f"Test frequencies: {frequencies}")
            print(f"nbins: {int(data['nbin'])}")

            freq, efstat = epoch_folding_search(
                np.array(EventList['time'].value),
                frequencies,
                nbin=int(data['nbin']),
                segment_size=float(data['segment_size']),
                gti=current_gti
            ) 
            # ---- PLOTTING --------
            plt.figure()
            plt.plot(freq, efstat, label='EF statistics', alpha=0.8)
            plt.axhline(int(data['nbin']) - 1, ls='--', lw=3, color='k', label='n - 1')
            plt.axvline(float(data['frequency']), lw=3, alpha=0.5, color='r', label='Correct frequency')
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('EF Statistics')
            _ = plt.legend()

            threshold = (np.max(efstat)*0.1+int(data['nbin']) - 1)
            best_x, best_y = search_best_peaks(freq,efstat,threshold)
            y_min, y_max = plt.ylim()
            offset = 0.04 * (y_max - y_min)

            '''
            for i, (x_value, y_value) in enumerate(zip(best_x, best_y)):
                label = 'peaks' if i == 0 else None  # Label only the first line
                plt.axvline(x_value, ls='dotted', lw=2, color='k', label=label)

                # Annotate the peak with its y-value slightly to the righ
                #plt.text(x_value + 0.03e-5, y_value, f"{y_value:.2f}", color='darkorange', fontsize=10)
                plt.text(x_value + 0.03e-5, y_value, f"{x_value:.4e}", color='darkorange', fontsize=10)
                plt.text(x_value + 0.03e-5, y_value-offset, f"{(float(data['frequency'])/x_value):.2f}", color='darkorange', fontsize=10)
            '''
            plt.savefig(output_plot)

            print(f"Plot saved: {output_plot}")
            print(f"Stingray version: {stingray.__version__}")
            print(f"epoch_folding_search is defined in: {inspect.getfile(epoch_folding_search)}")

            with h5py.File(output_file, 'w') as out:
                savehdf5(freq, efstat, out)

if __name__ == "__main__":
    main()
