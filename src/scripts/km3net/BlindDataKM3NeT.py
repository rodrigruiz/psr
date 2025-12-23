"""
Shuffle the location and energy information in HDF5 files to blind the data.

Usage: BlindDataKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files or file pattern
  -o --output_dir OUTPUT_DIR             Output directory

"""
# python3 BlindDataKM3NeT.py -i '/path/to/input/files/*.h5' -o /path/to/output/files/

from docopt import docopt
import os
import glob
from h5py import File
from astropy.table import Table
from astropy.time import Time
from astropy.io.misc.hdf5 import read_table_hdf5, write_table_hdf5
import numpy as np
from km3astro.io import load_hdf5_tables

def shuffle_columns(table, columns):
    """
    Shuffle the specified columns in the given Astropy table.

    Parameters
    ----------
    table : astropy.table.Table
        Table containing the columns to be shuffled.
    columns : list of str
        List of column names to shuffle.

    Returns
    -------
    astropy.table.Table
        New table with shuffled columns.
    """
    for column in columns:
        if column in table.colnames:
            data = table[column]
            shuffled_data = np.random.permutation(data)
            table[column] = shuffled_data
    return table

def main():
    arguments = docopt(__doc__)

    data = {key.replace("-", ""): value for key, value in arguments.items()}

    input_files = []
    for pattern in data['input_files']:
        input_files.extend(glob.glob(pattern))
    input_files.sort()

    if not input_files:
        print(f"No files matching pattern: {data['input_files']}")
        return

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    for file in input_files:
        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]
        output_filename = os.path.join(data['output_dir'], file_name + "_blinded.h5")

        tables = load_hdf5_tables(file)
        print(f"Processing file: {file}")

        # Check and shuffle times and columns for each table if it exists
        if hasattr(tables, 'id_table'):
             print("Shuffling times in ID table...")
             tables.id_table = shuffle_columns(tables.id_table, ['timeslice_utc_time'])

        #if hasattr(tables, 'mc_table') and tables.mc_table is not None:
        #    print("Shuffling columns in MC table...")
        #    tables.mc_table = shuffle_columns(tables.mc_table, ['energy', 'phi_detectorframe', 'theta_detectorframe'])

        #if hasattr(tables, 'reco_table'):
        #    print("Shuffling columns in RECO_EVENTS table...")
        #    tables.reco_table = shuffle_columns(tables.reco_table, ['energy', 'phi_detectorframe', 'theta_detectorframe'])

        # Write the shuffled tables back to an HDF5 file
        with File(output_filename, 'w') as h5file:
            # Write HEADER table
            if hasattr(tables, 'header_table'):
                write_table_hdf5(tables.header_table, h5file, path="HEADER", serialize_meta=True)
            # Write ID table
            if hasattr(tables, 'id_table'):
                write_table_hdf5(tables.id_table, h5file, path="ID", serialize_meta=True)
            # Write RECO_EVENTS table
            if hasattr(tables, 'reco_table'):
                reco_grp = h5file.create_group("RECO")
                write_table_hdf5(tables.reco_table, reco_grp, path="RECO_EVENTS", serialize_meta=True)
                # Write RECO_FITINF table
                write_table_hdf5(tables.fitinf_table, reco_grp, path="FITINF", serialize_meta=True)
            # Write MC_EVENTS table
            if hasattr(tables, 'mc_table') and tables.mc_table is not None:
                mc_grp = h5file.create_group("MC")
                write_table_hdf5(tables.mc_table, mc_grp, path="MC_EVENTS", serialize_meta=True)

        print(f"Shuffled file written: {output_filename}")

if __name__ == "__main__":
    main()