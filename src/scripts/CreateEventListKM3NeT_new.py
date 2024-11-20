""" Load KM3NeT root data files and convert them to astropy tables. 

Usage: CreateEventListKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR -s SOURCE_SPECS_FILE -e AR_SHOWER_FILE -t AR_TRACK_FILE [--energy_th=<float>] [--trackscore_th=<float>] [--detector=<detector>] [--energy_low=<energy_low>] [--energy_high=<energy_high>] [--shower_reco_name=<shower_reco_name>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files or file pattern
  -o --output_dir OUTPUT_DIR             Output directory
  -s --source SOURCE_SPECS_FILE          Hdf5 file containing information about the source of interest (ra, dec, P_orb, ...) 
  -e --ar_shower_file AR_SHOWER_FILE     Hdf5 file containing the angular resolution for shower like events
  -t --ar_track_file AR_TRACK_FILE       Hdf5 file containing the angular resolution for track like events
     --energy_th=<float>                 Energy Threshold. [default: 0]
     --trackscore_th=<float>             Track Score Threshold to choose track over shower reco. [default: 0.3]
     --detector=<string>                 Detector location ('arca','orca','antares') [default: arca]
     --energy_low=<int>                  Energy lower limit exponent (2 -> energy: 1e2) [default: 2]
     --energy_high=<int>                 Energy upper limit exponent (8 -> energy: 1e8) [default: 8]
     --shower_reco_name=<string>         Reco name of shower reco ('aashower' or 'jshower') [default: aashower]
"""
#python3 psr/src/scripts/CreateEventListKM3NeT.py -i '/home/hpc/capn/capn107h/software/hdf5TestOutput/*' -o eventlistTestOutput/ -s hdf5SourceFiles/Vela_X-1.h5

from docopt import docopt
import os, glob
import re

from h5py import File, Group
from astropy.table import Table, Column
from astropy.time import Time
from astropy.io.misc.hdf5 import read_table_hdf5
import astropy.units as u
from astropy.coordinates import SkyCoord
import km3io.definitions as kd

from km3astro.io import load_hdf5_tables

import numpy as np

from km3astro.coord import local_event
from km3astro import sources

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def extract_run_id(file_path):
    # Use regex to find the run_id after '.jterbr'
    match = re.search(r'\.jterbr(\d+)\.', file_path)
    
    if match:
        return match.group(1)
    else:
        return None


# Functions to fit to the angular resolution distribution
def logistic_function(E, a, b, c, d):
    return a / (1.0 + np.exp(-b * (np.log10(E) - c))) + d

def polynomial(E, a, b, c, d, e, f, g, h):
    return a * np.log10(E)**4 + b * np.log10(E)**3 + c * np.log10(E)**2 + d * np.log10(E) + e + f*np.log10(E)**5 +g * np.log10(E)**6 +h*np.log10(E)**7

# Function to read HDF5 file and fit the data to a model
def fit_angular_resolution(file_path, output_dir, energy_low, energy_high, fit_type="logistic"):
    # Open the HDF5 file
    table = Table.read(file_path, format='hdf5')

    # Access the columns by their names
    energy_bin_centers = table['energy_bin_center']
    angular_resolution = table['median_angular_separation']
    sigma1_width = table['sigma1_width']  # Assuming this is the uncertainty

    valid_indices = np.isfinite(angular_resolution) & np.isfinite(sigma1_width) & np.isfinite(energy_bin_centers)

    if not np.all(valid_indices):
        print("Warning: Some values were NaN or Inf and have been removed from the fit.")
        energy_bin_centers = energy_bin_centers[valid_indices]
        angular_resolution = angular_resolution[valid_indices]
        sigma1_width = sigma1_width[valid_indices]



    if fit_type == "logistic":
        # Initial guesses for the logistic function parameters
        popt, _ = curve_fit(logistic_function, energy_bin_centers, angular_resolution, 
                            p0=[3.5, 1, 4, 0.5], sigma=sigma1_width, absolute_sigma=True)
        fit_func = lambda E: logistic_function(E, *popt)
        plot_name = os.path.join(output_dir,f"fitted_angular_resolution_tracklike.png")
    elif fit_type == "polymial":
        # Initial guesses for the sinusoidal function parameters
        popt, _ = curve_fit(polynomial, energy_bin_centers, angular_resolution, sigma=sigma1_width, absolute_sigma=True)
        fit_func = lambda E: polynomial(E, *popt)
        plot_name = os.path.join(output_dir,f"fitted_angular_resolution_showerlike.png")
    else:
        raise ValueError("Unknown fit type. Choose 'logistic' or 'polymial'.")
    
    # Plot the data and the fit
    plt.figure(figsize=(10, 6))
    plt.plot(energy_bin_centers, angular_resolution, color='black', marker='o', linestyle='none', label='Data')
    energy_fit = np.logspace(energy_low, energy_high, 1000)  # Energy values for plotting the fit
    plt.plot(energy_fit, fit_func(energy_fit), color='firebrick', linestyle='-', linewidth=2, label='Fit')
    plt.xscale('log')
    # Add a shaded 1σ band for uncertainties (optional, depends on your data)
    sigma1_band_lower = angular_resolution - sigma1_width  # Example: replace sigma1_values with your actual 1σ uncertainties
    sigma1_band_upper = angular_resolution + sigma1_width
    plt.fill_between(energy_bin_centers, sigma1_band_lower, sigma1_band_upper, color='blue', alpha=0.3, label='1σ Band')

    # Add a shaded 2σ band for uncertainties (optional)
    sigma2_band_lower = angular_resolution - 2*sigma1_width   # Example: replace sigma2_values with your actual 2σ uncertainties
    sigma2_band_upper = angular_resolution + 2*sigma1_width
    plt.fill_between(energy_bin_centers, sigma2_band_lower, sigma2_band_upper, color='blue', alpha=0.1, label='2σ Band')


    plt.xscale('log')
    # Customize labels and title
    plt.xlabel('Energy [GeV]', fontsize=12)
    plt.ylabel('Angular Resolution [°]', fontsize=12)
    plt.title(f'Angular Resolution vs Energy ({fit_type} Fit Type)', fontsize=14)
    # Add gridlines for better readability
    plt.grid(True, which="both", linestyle="--", linewidth=0.7)
    plt.legend()
    plt.savefig(plot_name)

    return fit_func

def add_angular_resolution_to_events(event_table, shower_angres_function, track_angres_function, shower_reco_name = 'aashower'):
    """
    Calculate and add the angular resolution to the event_table based on the reconstruction type.

    Parameters:
    - event_table: The table containing event data.
    - shower_angres_function: Function to calculate angular resolution for shower events.
    - track_angres_function: Function to calculate angular resolution for track events.

    Returns:
    - event_table: Updated event_table with a new column for angular resolution.
    """
    
    # Create a new column for angular resolution with the same size as event_table
    angular_resolution = np.zeros(len(event_table))

    # Identify muon and shower event types
    muon_mask = (event_table['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & \
                (event_table['rec_stage'] == kd.reconstruction.JMUONBEGIN)
    
    if shower_reco_name == 'aashower':
        shower_mask = (event_table['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & \
                        (event_table['rec_stage'] == kd.reconstruction.AASHOWERBEGIN)
    elif shower_reco_name == 'jshower':
        shower_mask = (event_table['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & \
                    (event_table['rec_stage'] == kd.reconstruction.JSHOWERBEGIN)
    else: print("Unknwon shower reco type...")

    # Calculate angular resolutions using the vectorized functions
    angular_resolution[muon_mask] = track_angres_function(event_table['energy'][muon_mask])
    angular_resolution[shower_mask] = shower_angres_function(event_table['energy'][shower_mask])


    # Add the angular resolution column to the event_table
    # event_table['angular_resolution'] = angular_resolution

     # Add the angular resolution column to the event_table using add_column
    event_table.add_column(Column(name='angular_resolution', data=angular_resolution))


    return event_table

def select_events_based_on_track_score(tables, trackscore_threshold, shower_reco_name = "aashower" ):
    """
    Select events based on the track score and reconstruction type using a single mask.

    Parameters:
    - tables: The tables containing track score and reconstruction data.
    - trackscore_threshold: The threshold for track score to decide the selection.

    Returns:
    - event_table: A new table containing selected events based on the criteria.
    """

    print("Event Selection Track/Shower...")
    print(f"len(tables.id_table): {len(tables.id_table['event_id'])}")
    print(f"len(tables.reco_table): {len(tables.reco_table['event_id'])}")
    
    # Track score mask for high-scoring events
    high_score_mask = tables.id_table['track_score'] > trackscore_threshold

    # Create a dictionary mapping event_id to track_score mask (True for high score, False for low score)
    event_score_map = dict(zip(tables.id_table['event_id'], high_score_mask))

    # Create a mask that selects `muon` events for high score, `shower` events for low score
    combined_mask = np.zeros(len(tables.reco_table), dtype=bool)

    for i, event_id in enumerate(tables.reco_table['event_id']):
        if event_id in event_score_map:
            if event_score_map[event_id]:  # High score, select muon events
                if tables.reco_table['rec_type'][i] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE and \
                   tables.reco_table['rec_stage'][i] == kd.reconstruction.JMUONBEGIN:
                    combined_mask[i] = True
            else:  # Low score, select shower events
                if shower_reco_name == 'aashower':
                    if tables.reco_table['rec_type'][i] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE and \
                    tables.reco_table['rec_stage'][i] == kd.reconstruction.AASHOWERBEGIN:
                        combined_mask[i] = True
                elif shower_reco_name == 'jshower':
                    if tables.reco_table['rec_type'][i] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE and \
                    tables.reco_table['rec_stage'][i] == kd.reconstruction.JSHOWERBEGIN:
                        combined_mask[i] = True
                else: print("Wrong shower reco type specified (aashower or jshower implemented)... ")

    # Apply the combined mask to select the events
    selected_events = tables.reco_table[combined_mask]

    # Return the selected events as a new table
    event_table = Table(selected_events)

    return event_table


def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]

    input_files = []
    for pattern in data['input_files']:
        input_files.extend(glob.glob(pattern))
    # input_files.sort()

    if not input_files:
        print(f"No files matching pattern: {input_files}")
        return

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    energy_threshold = float(data['energy_th']) if data['energy_th'] != 0 else None
    trackscore_threshold = float(data['trackscore_th'])
    shower_reco_name = data['shower_reco_name']
    detector_location = data['detector']

    if detector_location == 'arca':
        shower_angres_function = fit_angular_resolution(data['ar_shower_file'], data['output_dir'], int(data['energy_low']), int(data['energy_high']), 'polymial')
        track_angres_function = fit_angular_resolution(data['ar_track_file'], data['output_dir'], int(data['energy_low']), int(data['energy_high']), 'logistic')
    elif detector_location == 'orca':
        shower_angres_function = fit_angular_resolution(data['ar_shower_file'], data['output_dir'], int(data['energy_low']), int(data['energy_high']), 'polymial')
        track_angres_function = fit_angular_resolution(data['ar_track_file'], data['output_dir'], int(data['energy_low']), int(data['energy_high']), 'polymial')
    else: print("Wrong detector name... ")

    with File(data['source'], 'r') as f:
        source_name = f['source_name'][()]
        print("Source Name:", source_name)

        # SkyCoord information
        ra = f['skycoord/ra'][()]
        dec = f['skycoord/dec'][()]
        skycoord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        print("SkyCoord:", skycoord)

    source_location = skycoord

    for file in input_files:
        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]

        run_id = extract_run_id(file_name)

        energy_threshold_str = str(energy_threshold) if energy_threshold is not None else "no"
        output_file = os.path.join(data['output_dir'], f"{file_name}_{energy_threshold_str}-energythreshold_eventlist_new")

        tables = load_hdf5_tables(file)

        event_table = select_events_based_on_track_score(tables, trackscore_threshold, shower_reco_name)

        print(event_table)

        event_table = add_angular_resolution_to_events(event_table, shower_angres_function, track_angres_function, shower_reco_name)

                # Find the events with angular resolution > 100 degrees
        problematic_events = event_table[event_table['angular_resolution'] > 100]

        print(f"Problematic Events (AngRes): {problematic_events['angular_resolution']}")
        print(f"Problematic Events (Energy): {problematic_events['energy']}")
        # Check the energy of these events
        # print(np.min(problematic_events['energy']))
        # print(np.max(problematic_events['energy']))
        # print(np.median(problematic_events['energy']))

        # Filter by reco_type
        # event_table = tables.reco_table[reco_type_mask]

        # Apply energy mask if necessary
        if energy_threshold is not None:
            energy_mask = event_table['energy'] >= energy_threshold
            event_table = event_table[energy_mask]

        # Extract the times
        # times = event_table['tracktime_utc'] if reco_type != "mc" else event_table['timeslice_utc_time']
        times = event_table['tracktime_utc']
        angular_resolution = event_table['angular_resolution']

        # Filter by source location (angular distance)
        if source_location is not None:
            event_location = local_event(np.array(event_table['theta_detectorframe']),
                                         np.array(event_table['phi_detectorframe']),
                                         times, detector_location)
            separation = event_location.separation(source_location)

            # print("Separation:", separation)
            # print("Angular resolution:", angular_resolution)


            

            location_mask = separation <= angular_resolution * u.deg

            # Apply location mask
            event_table = event_table[location_mask]

        # Ensure that after applying all the masks, the lengths of the columns are consistent
        # times = event_table['tracktime_utc'] if reco_type != "mc" else event_table['timeslice_utc_time']

        times = event_table['tracktime_utc']
        energy = event_table['energy']
        event_id = event_table['event_id']

        rec_type = event_table['rec_type']
        rec_stage = event_table['rec_stage']

        separation_deg = separation[location_mask].to(u.deg).value  # Filtered separation
        detector_name = [detector_location] * len(times)  # Detector name for each event
        run_ids = [run_id] * len(times)

        # Create the output table
        event_list = Table([ event_id, rec_type, rec_stage, times, energy, separation_deg, detector_name],
                           names=['event_id', 'rec_type', 'rec_stage', 'time', 'energy', 'separation', 'detector'])

        print(event_list)

        # Save the output table
        
        event_list.write(output_file + '.hdf5', format='hdf5', overwrite=True, serialize_meta=True)
        print("hdf5 File written...")



if __name__ == "__main__":
    main()