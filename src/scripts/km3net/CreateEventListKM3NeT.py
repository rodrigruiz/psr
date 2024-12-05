""" Load KM3NeT root data files and convert them to astropy tables. 

Usage: CreateEventListKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR -s SOURCE_SPECS_FILE [--reco_type=<reco_type>] [--format=<format>] [--energy_th=<float>] [--detector=<detector>] [--dist=<float>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files or file pattern
  -o --output_dir OUTPUT_DIR             Output directory
  -s --source SOURCE_SPECS_FILE          Hdf5 file containing information about the source of interest (ra, dec, P_orb, ...) 
     --reco_type=<string>                Reconstruction Type. [default: jmuon]
     --format=<string>                   Output Format. [default: hdf5]  
     --energy_th=<float>                 Energy Threshold. [default: 0]
     --detector=<string>                 Detector location ('arca','orca','antares') [default: arca]
     --dist=<float>                      Maximum angular source distance of included events [default: 5.0]
"""
#python3 psr/src/scripts/CreateEventListKM3NeT.py -i '/home/hpc/capn/capn107h/software/hdf5TestOutput/*' -o eventlistTestOutput/ -s hdf5SourceFiles/Vela_X-1.h5

from docopt import docopt
import os, glob

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

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    reco_type = data['reco_type']
    format = data['format']
    detector_location = data['detector']
    distance_to_source = float(data['dist'])

    energy_threshold = float(data['energy_th']) if data['energy_th'] != 0 else None

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
        energy_threshold_str = str(energy_threshold) if energy_threshold is not None else "no"
        output_file = os.path.join(data['output_dir'], f"{file_name}_{reco_type}-recotype_{energy_threshold_str}-energythreshold_eventlist")

        tables = load_hdf5_tables(file)

        # Apply the reco_type mask
        if reco_type == 'mc':
            event_table = tables.mc_table
        else:
            if reco_type == 'jmuon':
                reco_type_mask = (tables.reco_table['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & \
                                 (tables.reco_table['rec_stage'] == kd.reconstruction.JMUONBEGIN)
            elif reco_type == 'aashower':
                reco_type_mask = (tables.reco_table['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & \
                                 (tables.reco_table['rec_stage'] == kd.reconstruction.AASHOWERBEGIN)
            elif reco_type == 'jshower':
                reco_type_mask = (tables.reco_table['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & \
                                 (tables.reco_table['rec_stage'] == kd.reconstruction.JSHOWERBEGIN)
            else:
                raise ValueError(f"Unknown reco_type: {reco_type}")

            # Filter by reco_type
            event_table = tables.reco_table[reco_type_mask]

        # Apply energy mask if necessary
        if energy_threshold is not None:
            energy_mask = event_table['energy'] >= energy_threshold
            event_table = event_table[energy_mask]

        # Extract the times
        times = event_table['tracktime_utc'] if reco_type != "mc" else event_table['timeslice_utc_time']

        # Filter by source location (angular distance)
        if source_location is not None:
            event_location = local_event(np.array(event_table['theta_detectorframe']),
                                         np.array(event_table['phi_detectorframe']),
                                         times, detector_location)
            separation = event_location.separation(source_location)
            location_mask = separation <= distance_to_source * u.deg

            # Apply location mask
            event_table = event_table[location_mask]

        # Ensure that after applying all the masks, the lengths of the columns are consistent
        times = event_table['tracktime_utc'] if reco_type != "mc" else event_table['timeslice_utc_time']
        energy = event_table['energy']
        event_id = event_table['event_id']
        separation_deg = separation[location_mask].to(u.deg).value  # Filtered separation
        detector_name = [detector_location] * len(times)  # Detector name for each event

        # Create the output table
        event_list = Table([times, energy, event_id, separation_deg, detector_name],
                           names=['time', 'energy', 'event_id', 'separation', 'detector'])

        # Save the output table
        if format == 'hdf5':
            event_list.write(output_file + '.hdf5', format='hdf5', overwrite=True, serialize_meta=True)
        else:
            event_list.write(output_file + '.dat', format='ascii.ecsv', overwrite=True)


if __name__ == "__main__":
    main()