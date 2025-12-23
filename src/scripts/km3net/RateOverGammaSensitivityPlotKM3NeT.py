""" 

Usage: RateOverGammaSensitivityPlotKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR 

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

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from h5py import File
import pandas as pd

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
    output_file = os.path.join(output_dir, f"{common_prefix}_SensitvityRateOverGamma.hdf5")
    output_plot = os.path.join(output_dir, f"{common_prefix}_SensitvityRateOverGamma_plot.png")

    pval_dict = {}
    gamma_set = set()
    rate_set = set()
    gamma_5s_list = []
    r5sigma_list = []

    for file in input_files:
        with h5py.File(file, 'r') as h5_file:
            # Extract gamma
            p = Table.read(h5_file, path="sigma_thresholds")
            gamma = float(p['gamma'][0])
            r5sigma = float(p['SNR_5sigma'][0])
            gamma_set.add(gamma)
            gamma_5s_list.append(gamma)
            r5sigma_list.append(r5sigma)

            # Extract rates and p_values
            d = Table.read(h5_file, path="histogram_data")
            rates = np.array(d['SNR'])
            pvalues = np.array(d['p_Value'])

            # Populate dictionary
            for r, pv in zip(rates, pvalues):
                rate_set.add(r)
                pval_dict[(gamma, r)] = pv

    # Sort gamma and rate axes
    gamma_list = sorted(gamma_set)
    rate_list = sorted(rate_set)

    # Build 2D matrix
    pval_matrix = []
    for gamma in gamma_list:
        row = []
        for rate in rate_list:
            row.append(pval_dict.get((gamma, rate), None))
        pval_matrix.append(row)

    # Optional: convert to DataFrame
    df = pd.DataFrame(pval_matrix, index=gamma_list, columns=rate_list)
    df.index.name = "gamma"
    df.columns.name = "rate"


    
    # Save the combined EventList

    # Maybe save data:
    df.to_hdf(output_file, key='histogram_data', mode='w')


    # Transpose the p-value array so SNR is the vertical axis (rows)
    pval_array = df.to_numpy(dtype=float).T  # Transpose

    # Replace zero or tiny values to avoid log issues
    pval_array = np.clip(pval_array, 1e-20, None)

    # Axes: get index and columns
    gammas = df.index.to_numpy()   # X-axis
    snrs = df.columns.to_numpy()   # Y-axis

    # Create meshgrid for contours
    Gamma, SNR = np.meshgrid(gammas, snrs)

    # Create the plot
    plt.figure(figsize=(8, 6))
    c = plt.pcolormesh(gammas, snrs, pval_array, shading='auto', cmap='viridis', norm=LogNorm())

    #contour = plt.contour(Gamma, SNR, pval_array, levels=[2.87e-7], colors='red', linewidths=2)
    #plt.clabel(contour, fmt={2.87e-7: r'5$\sigma$ (2.87e-7)'}, inline=True, fontsize=10)7

        
    gamma_sorted, r5sigma_sorted = zip(*sorted(zip(gamma_5s_list, r5sigma_list)))

    # Plot the 5σ threshold line
    plt.plot(gamma_sorted, r5sigma_sorted, color='k', linestyle="--", linewidth=2, label=r'5$\sigma$ threshold')

    plt.colorbar(c, label='p-value (log scale)')
    plt.xlabel(r'Spectral Index $\gamma$')
    plt.ylabel(r'Signal Rate [1/d]')
    #plt.ylim([7,15])
    plt.title('p-value over signal rate and gamma')
    plt.tight_layout()
    plt.legend()
    plt.savefig(output_plot)

    print(f"Plot saved to {output_plot}")
    
    threshold_df = pd.DataFrame({'E_min': gamma_sorted, 'SNR_5sigma': r5sigma_sorted})
    threshold_df.to_hdf(output_file, key='sigma_thresholds', mode='a')

    
if __name__ == "__main__":
    main()