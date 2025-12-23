""" Fetching outputs from SignalNoiseStatisticsKM3NeT.py and creating a single data file.

Usage: CombineChi2ValuesKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR

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
    output_file = os.path.join(output_dir, f"{common_prefix}_ParameterStudy_5sigma_SNRs.hdf5")
    #output_plot = os.path.join(output_dir, f"{common_prefix}_SNR{ratio}_maxchi2_plot.png")

    max_chi2_list = []

    for file in input_files:
        with h5py.File(file, 'r') as h5_file:
            frequencies = h5_file['efstats/frequencies'][:]
            chi2s = h5_file['efstats/chi2s'][:]

            # Find the maximum chi-squared value and its corresponding frequency
            max_chi2_index = np.argmax(chi2s)
            max_chi2_value = chi2s[max_chi2_index]
            corresponding_frequency = frequencies[max_chi2_index]

            # Append the results to the list
            max_chi2_list.append((corresponding_frequency, max_chi2_value))



    
    max_chi2_values = [chi2 for _, chi2 in max_chi2_list]
    mean_max_chi2 = np.mean(max_chi2_values)
    meanvalue_data = np.array([ratio, mean_max_chi2])
    # Convert the list to an Astropy Table
    max_chi2_table = Table(rows=max_chi2_list, names=('Frequency', 'Max_Chi2'))
    #meanvalue_table = Table(rows=meanvalue_data, names=('SNR', 'Mean_Max_Chi2'))

    # Create a histogram of the maximum chi-squared values
    plt.figure()
    plt.hist(max_chi2_values, bins=nhbins, edgecolor='black')
    plt.axvline(mean_max_chi2, color='r', linestyle='dashed', linewidth=1)
    plt.title(f'Histogram of Maximum Chi-Squared Values (SNR={ratio})')
    plt.xlabel('Maximum Chi-Squared Value')
    plt.ylabel('Number of Entries')
    plt.grid(True)
    _ = plt.legend()
    plt.savefig(output_plot)   
    plt.close()     

    # Save the combined EventList
    max_chi2_table.write(output_file, format='hdf5', path='histogram_data', overwrite=True, serialize_meta = True)
    # meanvalue_table.write(output_file, format='hdf5', path='meanvalue_data', overwrite=True, serialize_meta = True)
    # Write the mean value to the HDF5 file using h5py
    with h5py.File(output_file, 'a') as h5_file:
        h5_file.create_dataset('meanvalue', data=meanvalue_data)

    print(f"Maximum Chi2 written to {output_file}")

    
if __name__ == "__main__":
    main()