"""
Classify KM3NeT events regarding track or shower-like signal.

Usage: 
    ClassifyEventsKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR

Options:
    -h --help                     Show this help message.
    -i --input_files INPUT_FILES   Input files or file pattern.
    -o --output_dir OUTPUT_DIR     Output directory.
"""

from docopt import docopt
import os, glob
import h5py
import numpy as np
from astropy.table import Column
from astropy.io.misc.hdf5 import write_table_hdf5
from km3astro.io import load_hdf5_tables


def calculate_track_score(id_table):
    """Generate a random track_score array based on the length of id_table."""
    return np.random.rand(len(id_table['event_id']))


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

    # Process each file
    for file in input_files:
        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]

        output_file_name = f"{file_name}_classified.h5"
        output_file_path = os.path.join(data['output_dir'], output_file_name)

        tables = load_hdf5_tables(file)
        print(f"Available tables in {file}: {list(tables.__dict__.keys())}")


        if hasattr(tables, 'id_table'):

            print("Calculating track score...")
            track_score = calculate_track_score(tables.id_table)
            track_score_column = Column(track_score, name= 'track_score')

            with h5py.File(output_file_path,'w') as h5file:
                
                tables.id_table.add_column(track_score_column)

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

