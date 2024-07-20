""" Fetching Chi2 distributions from Chi2HistogramKM3NeT.py and creating Chi2 over SNR plots.

Usage: SignalNoiseStatisticsKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR

Options:
  -h --help                              Show this help message
  -i --input_files INPUT_FILES...        Input files
  -o --output_dir OUTPUT_DIR             Output directory  
"""


from docopt import docopt
import os
import glob
import h5py
import numpy as np
from astropy.table import vstack, Table
import plens.EventList as EL
from matplotlib import pyplot as plt

def main():
    arguments = docopt(__doc__)

    input_files = arguments['--input_files']
    output_dir = arguments['--output_dir']

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Extract common prefix from input filenames
    common_prefix = os.path.commonprefix(input_files)
    # Remove any trailing non-alphanumeric characters from common prefix
    common_prefix = os.path.basename(common_prefix).rstrip("_-.")
    if not common_prefix:
        common_prefix = "Test"
    output_file = os.path.join(output_dir, f"{common_prefix}_StatisticOverSNR.hdf5")
    output_plot = os.path.join(output_dir, f"{common_prefix}_StatisticOverSNR_plot.png")

    max_chi2_list = []
    ratio_list = []

    for file in input_files:
        with h5py.File(file, 'r') as h5_file:
            mean_max_chi2 = h5_file['meanvalue'][1]
            ratio = h5_file['meanvalue'][0]

            # Append the results to the list
            max_chi2_list.append(mean_max_chi2)
            ratio_list.append(ratio)
    

    # Create a histogram of the maximum chi-squared values
    plt.figure()
    plt.scatter(ratio_list,max_chi2_list,lw=3)
    plt.title(f'Test Statistic')
    plt.xlabel('Signal to Noise Ratio (SNR)')
    plt.ylabel('Maximum Chi-Squared Value')
    plt.grid(True)
    _ = plt.legend()
    plt.savefig(output_plot)   
    plt.close()     

    # Convert the list to an Astropy Table
    chi2_over_snr_table = Table([ratio_list,max_chi2_list], names=('SNR', 'Max_Mean_Chi2'))

    # Save the combined EventList
    chi2_over_snr_table.write(output_file, format='hdf5', path='histogram_data', overwrite=True, serialize_meta = True)

    print(f"Maximum Chi2 written to {output_file}")

    
if __name__ == "__main__":
    main()
    