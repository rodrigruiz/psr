""" Fetching results from AngularResolutionKM3NeT.py and creating final anguar resolution over energy plots from multiple runs.

Usage: CombineAngularResolutionKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR [--runtype=<runtype>] [--recotype=<recotype>] [--pltscale=<pltscale>]

Options:
  -h --help                              Show this help message
  -i --input_files INPUT_FILES...        Input files
  -o --output_dir OUTPUT_DIR             Output directory  
     --runtype=<string>                  Run type ('nue','numu','anue','anumu') [default: anue]
     --recotype=<string>                 Reco type ('jmuon','aashower') [default: jmuon]
     --pltscale=<string>                 Scale of the y-axis ('linear','log') [default: linear]
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
    run_type = arguments['--runtype']
    reco_type = arguments['--recotype']
    plt_scale = arguments['--pltscale']

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Extract common prefix from input filenames
    common_prefix = os.path.commonprefix(input_files)
    # Remove any trailing non-alphanumeric characters from common prefix
    common_prefix = os.path.basename(common_prefix).rstrip("_-.")
    if not common_prefix:
        common_prefix = "Test"
    output_file = os.path.join(output_dir, f"{common_prefix}_AngResOverEnergy.hdf5")
    output_plot_lin = os.path.join(output_dir, f"{common_prefix}_AngResOverEnergy_plotlin.png")
    output_plot_log = os.path.join(output_dir, f"{common_prefix}_AngResOverEnergy_plotlog.png")
    
    energy_bin_center_ref = None
    combined_median = []
    combined_sigma1 = []
    counter = 0
    """
    for file in input_files:
        with h5py.File(file, 'r') as h5_file:
            energy_bin_center = h5_file['energy_bin_center']
            median_angular_separation = h5_file['median_angular_separation']

            if counter == 0: 
                energy_bin_center_list = energy_bin_center
                print("Energy bins initialy defined...")
            elif energy_bin_center != energy_bin_center_list:
                print("Error: Energy bins do NOT match with other files! They should...")
            
            median_angular_separation_list[counter] = median_angular_separation
    """
    
    for file in input_files:
        table = Table.read(file, format='hdf5')
        
            # Extract columns
        energy_bin_center = table['energy_bin_center']
        median_angular_separation = table['median_angular_separation']
        sigma1_width = table['sigma1_width']

        # Check if energy_bin_center is consistent
        if energy_bin_center_ref is None:
            energy_bin_center_ref = energy_bin_center
        elif not np.allclose(energy_bin_center_ref, energy_bin_center):
            raise ValueError(f"Energy bin centers do not match in file {file}")

        # Append the data from each file
        if counter == 0:
            # First file: Initialize the combined arrays
            combined_median = [median_angular_separation]
            combined_sigma1 = [sigma1_width]
        else:
            # Append subsequent files' data
            combined_median.append(median_angular_separation)
            combined_sigma1.append(sigma1_width)

    # Convert the combined lists to numpy arrays
    combined_median = np.array(combined_median)
    combined_sigma1 = np.array(combined_sigma1)  

    # Transpose
    combined_median = combined_median.T
    combined_sigma1 = combined_sigma1.T

    overall_median = np.mean(combined_median, axis=1)
    overall_sigma1 = np.mean(combined_sigma1, axis=1)

    print(energy_bin_center_ref)
    print(overall_median)
    print(overall_sigma1)

    bin_centers = np.array(energy_bin_center_ref)
    median_separations = overall_median
    sigma1_bands = overall_sigma1


    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(bin_centers, median_separations, color='black', label='Mean Angular Separation')
    plt.fill_between(bin_centers, median_separations - sigma1_bands, median_separations + sigma1_bands, color='blue', alpha=0.3, label='1σ Band')

    if plt_scale == 'linear': 
        plt.axhline(0,color="k",ls="dashed",lw=2)

    # plt.axhline(y=0.1, color='orange', linewidth=4)  # Example of a horizontal line
    plt.xscale('log')
    plt.yscale(plt_scale)
    #plt.ylim([1e-4,1e2])
    plt.xlabel('Energy [GeV]')
    plt.ylabel('Angular Resolution [°]')
    plt.title(f'Angular Resolution vs Energy for {reco_type} Reco Type ({run_type})')
    plt.legend()
    plt.grid(True, which="both", ls="--")

    output_plotname = "CombinedTestPlotAngularRes" + reco_type + "_" + run_type + "_" + plt_scale+ "_" + common_prefix + ".png"

    plt.savefig(output_plotname)   
    plt.close()   
    print(f"Plot saved: {output_plotname}")

    hdf5_name = "CombinedAngularResolutionOverEnergy_" +  reco_type + "_" + common_prefix + ".hdf5"
    angular_resolution_table = Table([bin_centers, median_separations, sigma1_bands], names=['energy_bin_center', 'median_angular_separation','sigma1_width'])
    angular_resolution_table.write(output_file, format='hdf5', overwrite=True, serialize_meta=True)
    print(f"hdf5 file saved: {hdf5_name}")
    
     # Plotting Histogram of All Separations
    plt.figure(figsize=(10, 6))
    plt.hist(overall_median, bins=53, color='darkorange')
    plt.yscale(plt_scale)
    plt.xlabel('Angular Separation [°]')
    plt.ylabel('Counts')
    plt.title(f'Histogram of Angular Separations for {reco_type} Reco Type ({run_type})')
    
    output_histname = "CombinedHistogramSeparations_" + reco_type + "_" + run_type + "_" + plt_scale + "_" + common_prefix + ".png"
    plt.savefig(output_histname)
    plt.close()
    print(f"Histogram saved: {output_histname}")
    

    
if __name__ == "__main__":
    main()
    