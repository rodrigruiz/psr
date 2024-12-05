"""
Classify KM3NeT events regarding track or shower-like signal.

Usage: 
    ClassifyEventsKM3NeT.py -i INPUT_FILES... -p=<parampid_folder> -o OUTPUT_DIR [--detector=<detector>]

Options:
    -h --help                                      Show this help message.
    -i --input_files INPUT_FILES                   Input files or file pattern.
    -p --parampid_folder=<parampid_folder>         Folder with the classified files from parampid
    -o --output_dir OUTPUT_DIR                     Output directory.
       --detector=<detector>                       Name of the Detector. 'arca' or 'orca' [default: arca]
"""

from docopt import docopt
import os, glob
import h5py
import numpy as np
from astropy.table import Column
from astropy.io.misc.hdf5 import write_table_hdf5
from km3astro.io import load_hdf5_tables
import pandas as pd


def calculate_random_track_score(id_table):
    """Generate a random track_score array based on the length of id_table."""
    return np.random.rand(len(id_table['event_id']))

def calculate_track_score(id_table, parampid_filepath):
    """Load track_score values from ParamPID file and map them to id_table based on event_id."""

    # Load the ParamPID file
    try:
        df = pd.read_hdf(parampid_filepath, 'summary')
    except Exception as e:
        raise IOError(f"Error reading ParamPID file: {parampid_filepath}. Details: {e}")

    # Ensure required columns exist
    if 'group_id' not in df.columns or 'track_score' not in df.columns:
        raise ValueError(f"ParamPID file {parampid_filepath} is missing required columns 'group_id' or 'track_score'.")


    print(f"Loaded ParamPID Data (first 5 rows):\n{df.head()}")

    # Map group_id (ParamPID) to event_id (id_table) for track_score
    event_ids = id_table['event_id']
    group_id_to_score = dict(zip(df['group_id'], df['track_score']))

    print(f"Mapping Dictionary (first 10 items): {list(group_id_to_score.items())[:10]}")
    print(f"Event IDs from ID Table (first 10): {event_ids[:10]}")
    
    # Find unmatched event_ids
    unmatched_event_ids = [event_id for event_id in event_ids if event_id not in group_id_to_score]

    if unmatched_event_ids:
        print(f"Warning: {len(unmatched_event_ids)} event_ids have no matching group_id in ParamPID file.")
        print(f"First 10 unmatched event_ids: {unmatched_event_ids[:10]}")
    else:
        print("All event_ids have matching group_id entries in the ParamPID file.")



    # Assign track scores to id_table, using NaN for missing group_ids
    track_scores = [group_id_to_score.get(event_id, np.nan) for event_id in event_ids]
    print(f"Mapped Track Scores (first 10): {track_scores[:10]}")

    return track_scores


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

    parampid_folder = str(data['parampid_folder'])
    det_name = str(data['detector'])

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    # Process each file
    for file in input_files:
        folder_path, file_name = os.path.split(file)

        end_of_name = f"_{det_name}.h5"

        if not file_name.endswith(end_of_name):
            raise ValueError(f"Input filename must end with '{end_of_name}'")


        file_name = os.path.splitext(file_name)[0]
        
        # base_name = file_name.replace(end_of_name, "")

        ## create filepaths for loading classified parampid files
        parampid_filename = file_name.replace(f"_{det_name}", ".root") + "_scored.h5"
        parampid_filepath = os.path.join(parampid_folder, parampid_filename)
        print(f"ParamPID Filepath: {parampid_filepath}")


        output_file_name = f"{file_name}_classified.h5"
        output_file_path = os.path.join(data['output_dir'], output_file_name)

        tables = load_hdf5_tables(file)
        print(f"Available tables in {file}: {list(tables.__dict__.keys())}")


        if hasattr(tables, 'id_table'):

            print("Calculating track score...")
            track_score = calculate_track_score(tables.id_table, parampid_filepath)
            track_score_column = Column(track_score, name= 'track_score')

            with h5py.File(output_file_path,'w') as h5file:
                
                tables.id_table.add_column(track_score_column)
                print(f"ID Table with Track Scores:\n{tables.id_table[:10]}")


                write_table_hdf5(tables.header_table, h5file, path='HEADER', serialize_meta=True)
                write_table_hdf5(tables.id_table, h5file, path='ID', serialize_meta=True)

                reco_grp = h5file.create_group("RECO")
                write_table_hdf5(tables.reco_table, reco_grp, path="RECO_EVENTS", serialize_meta=True)
                
                if tables.mc_table is not None:
                    mc_grp = h5file.create_group("MC")
                    write_table_hdf5(tables.mc_table, mc_grp, path='MC_EVENTS', serialize_meta=True)
                print(f"New HDF5 file written: {h5file}")
        else:
            print(f"Error loading 'id_table' from {file}.")


if __name__ == "__main__":
    main()

