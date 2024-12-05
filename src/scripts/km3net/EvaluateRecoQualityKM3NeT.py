""" Load KM3NeT hdf5 astropy tables and calculate the distribution of the angular resolution over the energy. 

Usage: AngularResolutionKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR [--detector=<detector>] [--runtype=<runtype>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files or file pattern
  -o --output_dir OUTPUT_DIR             Output directory
     --detector=<string>                 Detector location ('arca','orca','antares') [default: arca]
     --runtype=<string>                  Run type ('nue','numu','anue','anumu') [default: anue]
"""

from docopt import docopt
import os, glob
import matplotlib.pyplot as plt

from h5py import File, Group
from astropy.table import Table
from astropy.time import Time
from astropy.io.misc.hdf5 import read_table_hdf5
import astropy.units as u
from astropy.coordinates import SkyCoord
import km3io.definitions as kd
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d

from km3astro.io import load_hdf5_tables

import numpy as np

from km3astro.coord import local_event
from km3astro import sources


def make_lik_histogram(table, reco_type, reco_stage, reco_name, n_bins, run_type):
    # Apply mask for the specified reco_type and reco_stage
    mask = (table['rec_type'] == reco_type) & (table['rec_stage'] == reco_stage)
    filtered_table = table[mask]

    # Negate likelihood for jmuon reco type
    if reco_name == "jmuon":
        likelihoods = filtered_table['likelihood']
    else:
        likelihoods = -filtered_table['likelihood']

    hist, bin_edges = np.histogram(likelihoods, bins=n_bins, density=True)
    kde = gaussian_kde(likelihoods, bw_method='scott')

    x_grid = np.linspace(min(bin_edges), max(bin_edges), 1000)
    pdf = kde(x_grid)

    # Plotting Histogram of All Separations
    plt.figure(figsize=(10, 6))
    plt.hist(likelihoods, bins=n_bins, color='darkorange', label="Likelihood Distribution",density=True)
    plt.plot(x_grid, pdf, label='KDE Smoothed Distribution', color='black')
    plt.xlabel('Likelihood')
    plt.ylabel('Density')
    plt.title(f'Histogram of Likelihoods for {reco_name} Reco Type ({run_type}, {len(filtered_table)} events)')
    plt.legend()

    output_histname = "HistogramLikelihoods_" + reco_name + "_" + run_type + ".png"
    plt.savefig(output_histname)
    plt.close()

    return hist, bin_edges, kde


def precompute_cdf(kde, x_grid):
    cdf_values = np.array([kde.integrate_box_1d(-np.inf, x) for x in x_grid])
    return interp1d(x_grid, cdf_values, bounds_error=False, fill_value=(0.0, 1.0))

def complementary_cdf(likelihood, cdf_func):
    cdf_value = cdf_func(likelihood)
    return 1.0 - cdf_value

def favored_reco_type(table, kde_jm, kde_aas, x_grid_jm, x_grid_aas):
    favored_rows = {}

    cdf_jm_func = precompute_cdf(kde_jm, x_grid_jm)
    cdf_aas_func = precompute_cdf(kde_aas, x_grid_aas)

    # Iterate over each row in the reco_table
    for row in table:
        event_id = row['event_id']
        
        # Negate likelihood for jmuon reco type
        if row['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE and row['rec_stage'] == kd.reconstruction.JMUONBEGIN:
            likelihood = row['likelihood']
            P_jm = complementary_cdf(likelihood, cdf_jm_func)
            P_aas = 1.0
        elif row['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE and row['rec_stage'] == kd.reconstruction.AASHOWERBEGIN:
            likelihood = -row['likelihood']
            P_aas = complementary_cdf(likelihood, cdf_aas_func)
            P_jm = 1.0
        else:
            continue

        if event_id not in favored_rows:
            favored_rows[event_id] = row
        else:
            # Compare the current reco with the stored one
            stored_row = favored_rows[event_id]
            
            
            if stored_row['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE and stored_row['rec_stage'] == kd.reconstruction.JMUONBEGIN:
                stored_likelihood = stored_row['likelihood']
                stored_P_jm = complementary_cdf(stored_likelihood, cdf_jm_func)
                stored_P_aas = 1.0
            elif stored_row['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE and stored_row['rec_stage'] == kd.reconstruction.AASHOWERBEGIN:
                stored_likelihood = -stored_row['likelihood']
                stored_P_aas = complementary_cdf(stored_likelihood, cdf_aas_func)
                stored_P_jm = 1.0

            if (stored_row['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE and P_aas < stored_P_jm) or (stored_row['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE and P_jm < stored_P_aas):
                favored_rows[event_id] = row


    # Convert the favored_rows dictionary back to a table
    favored_table = Table(rows=favored_rows.values(), names=table.colnames)
    print(favored_table)

    return favored_table

def print_jmuon_aashower_ratio(favored_table):
    # Count the number of events where each reco type is favored
    jmuon_favored_count = np.sum((favored_table['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & (favored_table['rec_stage'] == kd.reconstruction.JMUONBEGIN))
    aashower_favored_count = np.sum((favored_table['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & (favored_table['rec_stage'] == kd.reconstruction.AASHOWERBEGIN))

    # Calculate the ratio
    ratio = jmuon_favored_count / aashower_favored_count if aashower_favored_count > 0 else np.inf

    # Print the ratio
    print(f"Ratio of jmuon to aashower favored reco types: {ratio:.2f}")
    print(f"jmuon favored: {jmuon_favored_count} events")
    print(f"aashower favored: {aashower_favored_count} events")

def main():
    arguments = docopt(__doc__)

    print(arguments)

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


    detector_location = data['detector']
    run_type = data['runtype']

    for file in input_files:

        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]

        output_file = data['output_dir'] + file_name + '_angularresolution.h5'
        print(output_file, os.path.exists(output_file))

        tables = load_hdf5_tables(file)

    n_bins = 53    

    hist_jm, bin_edges_jm, kde_jm= make_lik_histogram(tables.reco_table, reco_type=kd.reconstruction.JPP_RECONSTRUCTION_TYPE, reco_stage=kd.reconstruction.JMUONBEGIN, reco_name = "jmuon" , n_bins=n_bins,run_type=run_type)
    hist_aas, bin_edges_aas, kde_aas = make_lik_histogram(tables.reco_table, reco_type=kd.reconstruction.AANET_RECONSTRUCTION_TYPE, reco_stage=kd.reconstruction.AASHOWERBEGIN, reco_name = "aashower" , n_bins=n_bins,run_type=run_type)

    x_grid_jm = np.linspace(min(bin_edges_jm), max(bin_edges_jm), 1000)
    x_grid_aas = np.linspace(min(bin_edges_aas), max(bin_edges_aas), 1000)

    favored_table = favored_reco_type(tables.reco_table, kde_jm, kde_aas, x_grid_jm, x_grid_aas)
    print_jmuon_aashower_ratio(favored_table)
    favored_table.write(os.path.join(data['output_dir'], f"{file_name}_favored_reco.h5"), path='favored_reco_table', overwrite=True, serialize_meta=True)

if __name__ == "__main__":
    main()