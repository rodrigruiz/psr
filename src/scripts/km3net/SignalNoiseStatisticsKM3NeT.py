""" Fetching Chi2 distributions from Chi2HistogramKM3NeT.py and creating Chi2 over SNR plots.

Usage: SignalNoiseStatisticsKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR [--nbin=<nbin>]

Options:
  -h --help                              Show this help message
  -i --input_files INPUT_FILES...        Input files
  -o --output_dir OUTPUT_DIR             Output directory  
     --nbin=<int>                        Number of bins [default: 32]
"""


from docopt import docopt
import os
import glob
import h5py
import numpy as np
from astropy.table import vstack, Table
import plens.EventList as EL
from matplotlib import pyplot as plt
from scipy.stats.distributions import chi2

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
    output_plot_lin = os.path.join(output_dir, f"{common_prefix}_StatisticOverSNR_plotlin.png")
    output_plot_log = os.path.join(output_dir, f"{common_prefix}_StatisticOverSNR_plotlog.png")
    max_chi2_list = []
    ratio_list = []

    for file in input_files:
        with h5py.File(file, 'r') as h5_file:
            mean_max_chi2 = h5_file['meanvalue'][1]
            ratio = h5_file['meanvalue'][0]

            # Append the results to the list
            max_chi2_list.append(mean_max_chi2)
            ratio_list.append(ratio)

    #combined_list = list(zip(ratio_list, max_chi2_list))
    #combined_list.sort(key=lambda x: x[0])
    #sorted_ratio_list, sorted_max_chi2_list = zip(*combined_list)


    

    pvalue = chi2.sf(max_chi2_list , int(arguments['--nbin']) -1)

    # Create a histogram of the maximum chi-squared values
    # Plot with linear y-axis
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax2.scatter(ratio_list,max_chi2_list,alpha=0)
    ax1.scatter(ratio_list,pvalue,lw=3,label="Test Statistic")
    ax1.set_xlabel('Signal to Noise Ratio (SNR)')
    ax1.set_ylabel('p Value')
    ax2.set_ylabel(r'Maximum $\chi^2$')
    ax1.grid(True)
    ax1.axhline(y=2.87e-7, color='k', linestyle='--', linewidth=1, label=r"$5\sigma$ threshold")
    ax2.scatter(ratio_list,max_chi2_list,alpha=0)
    ax1.legend()
    # Save the linear plot
    plt.savefig(output_plot_lin)
    plt.close(fig)

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax2.scatter(ratio_list,max_chi2_list,alpha=0)
    ax1.scatter(ratio_list,pvalue,lw=3,label="Test Statistic")
    ax1.set_xlabel('Signal to Noise Ratio (SNR)')
    ax1.set_ylabel('p Value')
    ax2.set_ylabel(r'Maximum $\chi^2$')
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax1.grid(True)
    ax1.axhline(y=2.87e-7, color='k', linestyle='--', linewidth=1, label=f"$5\sigma$ threshold")
    ax2.scatter(ratio_list,max_chi2_list,alpha=0)
    ax1.legend()
    # Save the logarithmic plot
    plt.savefig(output_plot_log)
    plt.close(fig)




    # Convert the list to an Astropy Table
    chi2_over_snr_table = Table([ratio_list,max_chi2_list,pvalue], names=('SNR', 'Max_Mean_Chi2', 'p_Value'))

    # Save the combined EventList
    chi2_over_snr_table.write(output_file, format='hdf5', path='histogram_data', overwrite=True, serialize_meta = True)

    print(f"Maximum Chi2 written to {output_file}")

    
if __name__ == "__main__":
    main()
    