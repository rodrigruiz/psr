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

from km3astro.io import load_hdf5_tables

import numpy as np

from km3astro.coord import local_event
from km3astro import sources


def plot_angular_resolution_vs_energy(table, reco_type, reco_stage, reco_name, energy_range, n_bins, run_type):
    # Apply mask for the specified reco_type and reco_stage
    mask = (table['reco_type'] == reco_type) & (table['reco_stage'] == reco_stage)
    filtered_table = table[mask]

    print(len(filtered_table))
    
    # Logarithmic binning for energy
    energy_min, energy_max = energy_range
    bins = np.logspace(np.log10(energy_min), np.log10(energy_max), n_bins + 1)
    print(len(bins))
    
    # Prepare arrays to store results
    mean_separations = []
    sigma1_bands = []
    sigma2_bands = []
    bin_centers = []
    
    for i in range(len(bins) - 1):
        # Get events within the current energy bin
        bin_mask = (filtered_table['mc_energy'] >= bins[i]) & (filtered_table['mc_energy'] < bins[i+1])
        bin_data = filtered_table[bin_mask]['separation']
        
        if len(bin_data) > 0:
            # Calculate mean and standard deviations
            mean_sep = np.mean(bin_data)
            std_sep = np.std(bin_data)
            n = len(bin_data)
            mean_separations.append(mean_sep)
            sigma1_bands.append(std_sep / np.sqrt(n))  # 1-sigma uncertainty
            sigma2_bands.append(2 * std_sep / np.sqrt(n))  # 2-sigma uncertainty
            bin_centers.append(np.sqrt(bins[i] * bins[i+1]))  # Geometric mean for bin center
        else:
            mean_separations.append(np.nan)
            sigma1_bands.append(np.nan)
            sigma2_bands.append(np.nan)
            bin_centers.append(np.sqrt(bins[i] * bins[i+1]))
    
    # Convert to numpy arrays for easier plotting
    mean_separations = np.array(mean_separations)
    sigma1_bands = np.array(sigma1_bands)
    sigma2_bands = np.array(sigma2_bands)
    bin_centers = np.array(bin_centers)
    
    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(bin_centers, mean_separations, color='black', label='Mean Angular Separation')
    plt.fill_between(bin_centers, mean_separations - sigma1_bands, mean_separations + sigma1_bands, color='blue', alpha=0.3, label='1σ Band')
    plt.fill_between(bin_centers, mean_separations - sigma2_bands, mean_separations + sigma2_bands, color='blue', alpha=0.1, label='2σ Band')
    # plt.axhline(y=0.1, color='orange', linewidth=4)  # Example of a horizontal line
    plt.xscale('log')
    plt.ylim([0,100])
    plt.xlabel('Energy [GeV]')
    plt.ylabel('Angular Resolution [°]')
    plt.title(f'Angular Resolution vs Energy for {reco_name} Reco Type ({run_type}, {len(filtered_table)} events)')
    plt.legend()
    plt.grid(True, which="both", ls="--")

    output_plotname = "TestPlotAngularRes" + reco_name + "_" + run_type + ".png"

    plt.savefig(output_plotname)   
    plt.close()   

     # Plotting Histogram of All Separations
    plt.figure(figsize=(10, 6))
    plt.hist(filtered_table['separation'], bins=53, color='darkorange')
    plt.xlabel('Angular Separation [°]')
    plt.ylabel('Counts')
    plt.title(f'Histogram of Angular Separations for {reco_name} Reco Type ({run_type}, {len(filtered_table)} events)')
    
    output_histname = "HistogramSeparations_" + reco_name + "_" + run_type + ".png"
    plt.savefig(output_histname)
    plt.close()


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

        # print(tables.mc_table)
        # print(tables.reco_table)

        event_ids = []
        mc_energies = []
        reco_types = []
        reco_stages = []
        separations = []
        
        # Create a dictionary for quick lookup of MC events by event_id
        mc_dict = {row['event_id']: row for row in tables.mc_table}
        id_dict = {row['event_id']: row for row in tables.id_table}

        # Loop through each row in the RECO_EVENTS table
        for reco_event in tables.reco_table:
            event_id = reco_event['event_id']
            if event_id in mc_dict:
                mc_event = mc_dict[event_id]
                id_event = id_dict[event_id]
                
                # Calculate separation
                reco_theta = reco_event['theta_detectorframe']
                reco_phi = reco_event['phi_detectorframe']
                mc_theta = mc_event['theta_detectorframe']
                mc_phi = mc_event['phi_detectorframe']

                mc_location = local_event(np.array(mc_event['theta_detectorframe']),np.array(mc_event['phi_detectorframe']),id_event['timeslice_utc_time'],detector_location)
                reco_location = local_event(np.array(reco_event['theta_detectorframe']),np.array(reco_event['phi_detectorframe']),reco_event['tracktime_utc'],detector_location)
                separation = reco_location.separation(mc_location)
                

                #mc_object = SkyCoord(np.array(mc_event['theta_detectorframe'])*u.rad,np.array(mc_event['phi_detectorframe'])*u.rad)
                #reco_object = SkyCoord(np.array(reco_event['theta_detectorframe'])*u.rad,np.array(reco_event['phi_detectorframe'])*u.rad)
                #separation = reco_object.separation(mc_object)
                #print(separation)
            
                # Store results
                event_ids.append(event_id)
                mc_energies.append(mc_event['energy'])
                reco_types.append(reco_event['rec_type'])
                reco_stages.append(reco_event['rec_stage'])
                separations.append(separation)

        result_table = Table([event_ids, mc_energies*u.GeV, reco_types, reco_stages, separations*u.deg],
                         names=('event_id', 'mc_energy', 'reco_type', 'reco_stage', 'separation'))
        print(result_table)

    energy_range = (1e3, 1e7)  # [GeV]
    n_bins = 20
    plot_angular_resolution_vs_energy(result_table, reco_type=101, reco_stage=300, reco_name = "jmuon" , energy_range=energy_range, n_bins=n_bins,run_type=run_type)
    plot_angular_resolution_vs_energy(result_table, reco_type=4000, reco_stage=0, reco_name = "aashower" , energy_range=energy_range, n_bins=n_bins,run_type=run_type)

if __name__ == "__main__":
    main()