""" Load KM3NeT root data files and convert them to astropy tables. 

Usage: CreateEventListKM3NeT_new.py -i INPUT_FILES... -o OUTPUT_DIR -s SOURCE_SPECS_FILE -e AR_SHOWER_FILE -t AR_TRACK_FILE [--energy_th=<float>] [--trackscore_th=<float>] [--muonscore_th=<float>] [--detector=<detector>] [--energy_low=<energy_low>] [--energy_high=<energy_high>] [--delta_search_min=<delta_search_min>] [--ang_res_source=<ang_res_source>] [--cone_all=<cone_all>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files or file pattern
  -o --output_dir OUTPUT_DIR             Output directory
  -s --source SOURCE_SPECS_FILE          Hdf5 file containing information about the source of interest (ra, dec, P_orb, ...) 
  -e --ar_shower_file AR_SHOWER_FILE     Hdf5 file containing the angular resolution for shower like events
  -t --ar_track_file AR_TRACK_FILE       Hdf5 file containing the angular resolution for track like events
     --energy_th=<float>                 Energy Threshold. [default: 0]
     --trackscore_th=<float>             Track Score Threshold to choose track over shower reco. [default: 0.3]
     --muonscore_th=<float>              Muon Score Threshold to flag events as muons and ignore them. [default: 0.9]
     --detector=<string>                 Detector location ('arca','orca','antares') [default: arca]
     --energy_low=<int>                  Energy lower limit exponent (2 -> energy: 1e2) [default: 2]
     --energy_high=<int>                 Energy upper limit exponent (8 -> energy: 1e8) [default: 8]
     --delta_search_min=<float>          Minimal angular search cone size in degrees [default: 8]
     --ang_res_source=<float>            Angular resolution of source in marcsec [default: 1]
     --cone_all=<cone_all>               Whether or not to include all events inside the maximum selection angle, independent of the angular resolution ('True' or 'False') [default: False] 
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

from km3astro.coord import local_event, neutrino_to_source_direction
from km3astro import sources

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

    
def extract_run_id(filename: str) -> str:
    # Find all 8-digit numbers in the filename
    matches = re.findall(r'(?<!\d)\d{8}(?!\d)', filename)
    
    # Return the last occurrence if there are multiple
    if matches:
        return matches[-1]  # Last occurrence
    
    return None  # No 8-digit number found


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

def add_angular_resolution_to_events(event_table, shower_angres_function, track_angres_function, detector_name = 'arca'):
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
    
    if detector_name == 'arca':
        shower_mask = (event_table['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & \
                        (event_table['rec_stage'] == kd.reconstruction.AASHOWERBEGIN)
    elif detector_name == 'orca':
        shower_mask = (event_table['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & \
                    (event_table['rec_stage'] == kd.reconstruction.JSHOWERBEGIN)
    else: print("Unknwon shower reco type...")

    # Calculate angular resolutions using the vectorized functions
    angular_resolution[muon_mask] = track_angres_function(event_table['energy'][muon_mask])
    angular_resolution[shower_mask] = shower_angres_function(event_table['energy'][shower_mask])


    # Add the angular resolution column to the event_tabl
    # event_table['angular_resolution'] = angular_resolution

     # Add the angular resolution column to the event_table using add_column
    event_table.add_column(Column(name='angular_resolution', data=angular_resolution*u.deg))


    return event_table

def select_events_based_on_track_score_old(tables, trackscore_threshold, detector_name = "arca" ):
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
                if detector_name == 'arca':
                    if tables.reco_table['rec_type'][i] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE and \
                    tables.reco_table['rec_stage'][i] == kd.reconstruction.AASHOWERBEGIN:
                        combined_mask[i] = True
                elif detector_name == 'orca':
                    if tables.reco_table['rec_type'][i] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE and \
                    tables.reco_table['rec_stage'][i] == kd.reconstruction.JSHOWERBEGIN:
                        combined_mask[i] = True
                else: print("Wrong shower reco type specified (aashower or jshower implemented)... ")

    # Apply the combined mask to select the events
    selected_events = tables.reco_table[combined_mask]

    # Return the selected events as a new table
    event_table = Table(selected_events)

    return event_table

def select_events_old(tables, trackscore_threshold, muonscore_threshold, detector_name="arca", is_mc = False):
    """
    Select events based on the track score and reconstruction type using a single mask.
    Also stores the track_score and muon_score in the selected event table.

    Parameters:
    - tables: The tables containing track score and reconstruction data.
    - trackscore_threshold: The threshold for track score to decide the selection.
    - detector_name: The detector name (arca or orca).

    Returns:
    - event_table: A new table containing selected events with track_score and muon_score.
    """

    print("Event Selection Track/Shower...")
    print(f"len(tables.id_table): {len(tables.id_table['event_id'])}")
    print(f"len(tables.reco_table): {len(tables.reco_table['event_id'])}")
    
    # Create a dictionary mapping event_id to (track_score, muon_score)
    #event_scores = {eid: (ts, ms) for eid, ts, ms in zip(tables.id_table['event_id'], 
    #                                                     tables.id_table['track_score'], 
    #                                                     tables.id_table['muon_score'])}
    event_scores = {eid: ts for eid, ts in zip(tables.id_table['event_id'], 
                                                         tables.id_table['track_score'])}
    # Track score mask for high-scoring events
    high_score_mask = tables.id_table['track_score'] > trackscore_threshold
    event_score_map = dict(zip(tables.id_table['event_id'], high_score_mask))

    # Create a mask that selects muon events for high score, shower events for low score
    combined_mask = np.zeros(len(tables.reco_table), dtype=bool)
    selected_event_ids = []

    for i, event_id in enumerate(tables.reco_table['event_id']):
        if event_id in event_score_map:
            if event_score_map[event_id]:  # High score, select muon events
                if tables.reco_table['rec_type'][i] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE and \
                   tables.reco_table['rec_stage'][i] == kd.reconstruction.JMUONBEGIN:
                    combined_mask[i] = True
                    selected_event_ids.append(event_id)
            else:  # Low score, select shower events
                if detector_name == 'arca':
                    if tables.reco_table['rec_type'][i] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE and \
                       tables.reco_table['rec_stage'][i] == kd.reconstruction.AASHOWERBEGIN:
                        combined_mask[i] = True
                        selected_event_ids.append(event_id)
                elif detector_name == 'orca':
                    if tables.reco_table['rec_type'][i] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE and \
                       tables.reco_table['rec_stage'][i] == kd.reconstruction.JSHOWERBEGIN:
                        combined_mask[i] = True
                        selected_event_ids.append(event_id)
                else:
                    print("Wrong shower reco type specified (aashower or jshower implemented)... ")

    # Apply the combined mask to select the events
    selected_events = tables.reco_table[combined_mask]

    # Hits threshold maybe her like this
    n_hits_threshold = 1
    n_hits_mask = tables.fitinf_table['initial_number_of_hits'] > n_hits_threshold

    print("tables.fitinf_table['initial_number_of_hits']: ",tables.fitinf_table['initial_number_of_hits'])
    print("n_hits_mask: ", n_hits_mask)
    print("len(n_hits_mask): ", len(n_hits_mask))

    hits_event_map = dict(zip(tables.fitinf_table['event_id'], n_hits_mask))
    print("hits_event_map: ", hits_event_map)

    #hitsthreshold_selected_events = []
    #for row in selected_events:
    #    event_id = row['event_id']
    #   if hits_event_map.get(event_id, False):  # Only include if event_id is in map and passes threshold
    #        hitsthreshold_selected_events.append(row)

    mask = np.array([hits_event_map.get(event_id, False) for event_id in selected_events['event_id']])

    # Step 2: Apply it directly to get the final selection
    hitsthreshold_selected_events = selected_events[mask]

    print("hitsthreshold_selected_Events: ", hitsthreshold_selected_events)

    # Retrieve corresponding track_score and muon_score
    track_scores = [event_scores[eid] for eid in selected_event_ids]
    #muon_scores = [event_scores[eid][1] for eid in selected_event_ids]

    # Create event table and add new columns
    event_table = Table(hitsthreshold_selected_events)
    event_table['track_score'] = track_scores
    #event_table['muon_score'] = muon_scores

    if is_mc:
        mc_weights_map = dict(zip(tables.mc_table['event_id'], tables.mc_table['normalized_weight']))
        event_weights = [mc_weights_map.get(eid, None) for eid in event_table['event_id']]  # Default weight 1.0 if missing
        event_table['normalized_weight'] = event_weights

    #muon_score_mask = event_table['muon_score'] <= muonscore_threshold
    #event_table = event_table[muon_score_mask]

    return event_table

def select_events(tables, trackscore_threshold, muonscore_threshold, detector_name="arca", is_mc = False):
    """
    Select events based on the track score and reconstruction type using a single mask.
    Also stores the track_score and muon_score in the selected event table.

    Parameters:
    - tables: Namespace-like object with .id_table, .reco_table, .fitinf_table, and optionally .mc_table
    - trackscore_threshold: Threshold above which the event is considered a track (muon), below is shower.
    - muonscore_threshold: (Currently unused) Placeholder for applying further muon score filtering.
    - detector_name: Either "arca" or "orca", to determine which shower reco type to expect.
    - is_mc: If True, add normalized MC weights to the final output.

    Returns:
    - event_table: A table with selected events, including track_score and all relevant fitinf columns.
    """

    print("Event Selection: Track/Shower Classification...")

    # Create masks for track_score above threshold
    id_event_ids = tables.id_table['event_id']
    track_scores = tables.id_table['track_score']
    high_score_mask = track_scores > trackscore_threshold
    event_score_map = dict(zip(id_event_ids, high_score_mask))
    score_lookup = dict(zip(id_event_ids, track_scores))  # For later column assignment

    # Loop through reco_table and select events by rec_type/stage based on track_score classification
    selected_event_ids = []
    reco_event_ids = tables.reco_table['event_id']
    combined_mask = np.zeros(len(reco_event_ids), dtype=bool)

    for i, event_id in enumerate(reco_event_ids):
        if event_id not in event_score_map:
            continue
        
        is_track = event_score_map[event_id]
        rec_type = tables.reco_table['rec_type'][i]
        rec_stage = tables.reco_table['rec_stage'][i]

        if is_track:
            # Select muon (track-like) events
            if rec_type == kd.reconstruction.JPP_RECONSTRUCTION_TYPE and rec_stage == kd.reconstruction.JMUONBEGIN:
                combined_mask[i] = True
                selected_event_ids.append(event_id)
        else:
            # Select shower-like events
            if detector_name == 'arca' and rec_type == kd.reconstruction.AANET_RECONSTRUCTION_TYPE and rec_stage == kd.reconstruction.AASHOWERBEGIN:
                combined_mask[i] = True
                selected_event_ids.append(event_id)
            elif detector_name == 'orca' and rec_type == kd.reconstruction.JPP_RECONSTRUCTION_TYPE and rec_stage == kd.reconstruction.JSHOWERBEGIN:
                combined_mask[i] = True
                selected_event_ids.append(event_id)

    selected_events = tables.reco_table[combined_mask]

    # Apply hit threshold using fitinf_table
    hit_mask = tables.fitinf_table['initial_number_of_hits'] > 1
    fitinf_map = {eid: row for eid, row in zip(tables.fitinf_table['event_id'], tables.fitinf_table[hit_mask])}

    # Filter selected events by hit threshold
    final_mask = np.array([event_id in fitinf_map for event_id in selected_events['event_id']])
    hit_filtered_events = selected_events[final_mask]

    # Reconstruct final event table
    event_table = Table(hit_filtered_events)

    # Add track_score column
    event_table['track_score'] = [score_lookup[eid] for eid in event_table['event_id']]

    # Add all relevant fitinf columns
    fitinf_columns = tables.fitinf_table.colnames
    for col in fitinf_columns:
        if col == 'event_id':
            continue
        event_table[col] = [fitinf_map[eid][col] for eid in event_table['event_id']]

    # Add MC weights if applicable
    if is_mc:
        weight_map = dict(zip(tables.mc_table['event_id'], tables.mc_table['normalized_weight']))
        event_table['normalized_weight'] = [weight_map.get(eid, 1.0) for eid in event_table['event_id']]

    return event_table


def calc_search_cone(ang_res_km3net, min_err, ang_res_source, cone_all):
    # Ensure all inputs have compatible units
    ang_res_km3net = ang_res_km3net.to(u.deg)
    min_err = min_err.to(u.deg)
    ang_res_source = ang_res_source.to(u.deg)
    
    # Extract the numerical values for computation
    ang_res_km3net_val = ang_res_km3net.value
    min_err_val = min_err.value
    ang_res_source_val = ang_res_source.value
    
    # Perform the calculation
    if cone_all is True:
        search_cone_val = np.pi/2 * np.sqrt(min_err_val**2 + ang_res_source_val**2)
    else:
        search_cone_val = np.pi/2 * np.sqrt(np.minimum(ang_res_km3net_val**2, min_err_val**2) + ang_res_source_val**2) 
    
    # Reapply the units (deg) to the result
    return search_cone_val * u.deg


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
    muonscore_threshold  = float(data['muonscore_th'])
    detector_name = data['detector']
    min_err = float(data['delta_search_min'])/1.58 * u.deg
    ang_res_source = float(data['ang_res_source']) * u.marcsec
    cone_all = str(arguments['--cone_all'])
    cone_all_bool = cone_all == 'True'


    if detector_name == 'arca':
        shower_angres_function = fit_angular_resolution(data['ar_shower_file'], data['output_dir'], int(data['energy_low']), int(data['energy_high']), 'polymial')
        track_angres_function = fit_angular_resolution(data['ar_track_file'], data['output_dir'], int(data['energy_low']), int(data['energy_high']), 'logistic')
    elif detector_name == 'orca':
        shower_angres_function = fit_angular_resolution(data['ar_shower_file'], data['output_dir'], int(data['energy_low']), int(data['energy_high']), 'polymial')
        track_angres_function = fit_angular_resolution(data['ar_track_file'], data['output_dir'], int(data['energy_low']), int(data['energy_high']), 'logistic')
    else: print("Wrong detector name... ")

    with File(data['source'], 'r') as f:
        source_name = f['source_name'][()].decode('utf-8')
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
        

        tables = load_hdf5_tables(file)
        if hasattr(tables,"mc_table") and tables.mc_table is not None:
            is_mc = True
        else: is_mc = False

        event_table = select_events(tables, trackscore_threshold, muonscore_threshold, detector_name, is_mc)

        #print(event_table)

        event_table = add_angular_resolution_to_events(event_table, shower_angres_function, track_angres_function, detector_name)

                # Find the events with angular resolution > 100 degrees
        problematic_events = event_table[event_table['angular_resolution'] > 45]
        print(f"Number of Events with Agnular Resolution > 45: {len(problematic_events['angular_resolution'])}")
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

        nanmask = np.isfinite(event_table['theta_detectorframe']) & np.isfinite(event_table['phi_detectorframe'])  # Add more columns if necessary
        event_table = event_table[nanmask]

        # Extract the times
        # times = event_table['tracktime_utc'] if reco_type != "mc" else event_table['timeslice_utc_time']
        #times = event_table['tracktime_utc']
        #angular_resolution = event_table['angular_resolution']

        search_cone_val = calc_search_cone(event_table['angular_resolution'], min_err, ang_res_source, cone_all_bool )
        # Mask rows with NaN values in any column
        

        # Filter by source location (angular distance)
        if source_location is not None:

            theta_array = np.array(event_table['theta_detectorframe'])
            phi_array = np.array(event_table['phi_detectorframe'])

            #print("Theta: ", theta_array)
            #print("Phi: ", phi_array)

            if theta_array.size > 0:
                print("Theta Min: ", np.min(theta_array))
                print("Theta Max: ", np.max(theta_array))
            else:
                print("Error: Theta array is empty")

            if phi_array.size > 0:
                print("Phi Min: ", np.min(phi_array))
                print("Phi Max: ", np.max(phi_array))
            else:
                print("Error: Phi array is empty")

            event_location = local_event(np.array(event_table['theta_detectorframe']),
                                         np.array(event_table['phi_detectorframe']),
                                         event_table['tracktime_utc'], detector_name)
            separation = event_location.separation(source_location)

            print("separation: ", separation)
            print("separation min: ", np.min(separation))
            print("separation max: ", np.max(separation))

            location_mask = separation <= search_cone_val # angular_resolution #*u.deg 
            # (if it doesn't work, remove the *u.deg in the creation of the angular resolution Table above and include it here)

            # Apply location mask
            event_table = event_table[location_mask]

        # Ensure that after applying all the masks, the lengths of the columns are consistent
        # times = event_table['tracktime_utc'] if reco_type != "mc" else event_table['timeslice_utc_time']

        times = event_table['tracktime_utc']
        energy = event_table['energy']
        event_id = event_table['event_id']
        track_score = event_table['track_score']
        #muon_score = event_table['muon_score']

        if is_mc is True: normalized_weight = event_table['normalized_weight']
                                


        rec_type = event_table['rec_type']
        rec_stage = event_table['rec_stage']

        print("##############################")
        print(f"Amount of Events: {len(times)}")

        #theta = event_table['theta_detectorframe']
        #phi = event_table['phi_detectorframe']
        theta = np.asarray(event_table['theta_detectorframe'], dtype=np.float64)
        phi = np.asarray(event_table['phi_detectorframe'], dtype=np.float64)

        if np.any(np.isnan(theta)):
            print("There are NaN values in theta!")
        else:
            print("No NaN values in theta.")

        if np.any(np.isnan(phi)):
            print("There are NaN values in phi!")
        else:
            print("No NaN values in phi.")

        print("Theta Max:", np.max(theta))
        print("Theta Min:", np.min(theta))
        print("Theta dtype:", theta.dtype)


        print(np.pi)
        print(np.all(theta <= np.pi))
        print(phi)
        print("Theta values exceeding π:", theta[theta > np.pi])
        assert np.all(theta <= np.pi), f"Error: Some theta values exceed π! Max: {np.max(theta)}"

        #theta = np.asarray(theta, dtype=np.float64)
        #phi = np.asarray(phi, dtype=np.float64)

        #source_azimuths , source_zeniths = neutrino_to_source_direction(np.array(theta*u.rad), np.array(phi*u.rad),radian=True)
        source_azimuths, source_zeniths = neutrino_to_source_direction(phi, theta, radian=True)


        separation_deg = separation[location_mask].to(u.deg).value  # Filtered separation
        detector_name_list = [detector_name] * len(event_table['tracktime_utc'])  # Detector name for each event
        run_ids = [run_id] * len(event_table['tracktime_utc'])

        theta = np.asarray(event_table['theta_detectorframe'], dtype=np.float64)
        phi = np.asarray(event_table['phi_detectorframe'], dtype=np.float64)

        '''

        # Create the output table
        if is_mc:
            event_list = Table([run_ids, event_id , rec_type, rec_stage, times, energy, theta , phi, source_azimuths, source_zeniths, separation_deg, track_score, normalized_weight, detector_name_list],
                           names=['run_id','event_id', 'rec_type', 'rec_stage', 'time','energy', 'theta_detectorframe', 'phi_detectorframe', 'azimuth', 'zenith', 'separation','track_score', 'normalized_weight', 'detector'],
                           dtype=['int64', 'int64', 'int64', 'int64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'str'] )
        else:
            event_list = Table([run_ids, event_id , rec_type, rec_stage, times, energy, theta , phi, source_azimuths, source_zeniths, separation_deg, track_score, detector_name_list],
                                names=['run_id','event_id', 'rec_type', 'rec_stage', 'time','energy', 'theta_detectorframe', 'phi_detectorframe', 'azimuth', 'zenith', 'separation','track_score', 'detector'],
                                dtype=['int64', 'int64', 'int64', 'int64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'str'] )
        '''
        # Add normalized_weight only if it's MC and the column exists
        if is_mc and 'normalized_weight' not in event_table.colnames:
            print("Warning: is_mc=True but 'normalized_weight' not in event_table columns!")
        elif not is_mc and 'normalized_weight' in event_table.colnames:
            # Drop the column for real data to avoid inconsistencies
            event_table.remove_column('normalized_weight')

        # Rename 'tracktime_utc' to 'time'
        
        event_list = event_table.copy()
        if 'tracktime_utc' in event_list.colnames:
            event_list.rename_column('tracktime_utc', 'time')
        else:
            raise KeyError("'tracktime_utc' column not found in event_table")
        
        # Add extra columns
        event_list['run_id'] = run_ids
        event_list['detector'] = detector_name_list
        event_list['azimuth'] = source_azimuths
        event_list['zenith'] = source_zeniths
        event_list['separation'] = separation_deg

        print(event_list)
        min_opening_angle = str(data['delta_search_min'])

        output_file = os.path.join(data['output_dir'], f"{file_name}_{source_name}_{min_opening_angle}_deg_{energy_threshold_str}-E_th_eventlist_new")
        # Save the output table
        
        
        event_list.write(output_file + '.hdf5', format='hdf5', overwrite=True, serialize_meta=True)
        print("hdf5 File written...")



if __name__ == "__main__":
    main()