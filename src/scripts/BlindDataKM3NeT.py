"""
Shuffle the timing information in HDF5 files to blind the data.

Usage: BlindDataKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files or file pattern
  -o --output_dir OUTPUT_DIR             Output directory

"""
#python3 BlindDataKM3NeT.py -i '/path/to/input/files/*.h5' -o /path/to/output/files/


from docopt import docopt
import os, glob

from h5py import File, Group
from astropy.table import Table
from astropy.time import Time
from astropy.io.misc.hdf5 import read_table_hdf5, write_table_hdf5
import astropy.units as u
from astropy.coordinates import SkyCoord
import km3io.definitions as kd

from km3astro.io import load_hdf5_tables

import numpy as np

from km3astro.coord import local_event
from km3astro import sources

def shuffle_times(table, time_column):
    """
    Shuffle the times in the given Astropy table.

    Parameters
    ----------
    table : astropy.table.Table
        Table containing the times to be shuffled.
    time_column : str
        The name of the time column to shuffle.

    Returns
    -------
    astropy.table.Table
        New table with shuffled times.
    """
    times = table[time_column]
    shuffled_times = np.random.permutation(times)
    table[time_column] = shuffled_times
    return table

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
    

    for file in input_files:
        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]
        output_filename = os.path.join(data['output_dir'], file_name + "_blinded.h5")

        tables = load_hdf5_tables(file)
        #print(tables.__dict__)
        #print(tables)
        shuffled_tables = {}

        # Check and shuffle times for each reco type table if it exists
        reco_types = {
            'mc': ('id_table', 'timeslice_utc_time'),
            'reco': ('reco_table', 'tracktime_utc'),
        }

        for reco_type, details in reco_types.items():
            print(reco_type)
            print(details)
            if reco_type == 'mc' and hasattr(tables, 'mc_table') and tables.mc_table is not None:
                print("MC times are present and will now get shuffled...")
                table_name, time_column = details
                table = tables.id_table
                shuffled_tables[table_name] = shuffle_times(table, time_column)
            elif reco_type == 'reco':
                print("RECO times are present and will now get shuffled...")
                table_name, time_column = details
                table = tables.reco_table
                shuffled_tables[table_name] = shuffle_times(table, time_column)
        
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
            # Write MC_EVENTS table
            if hasattr(tables, 'mc_table') and tables.mc_table is not None:
                mc_grp = h5file.create_group("MC")
                write_table_hdf5(tables.mc_table, mc_grp, path="MC_EVENTS", serialize_meta=True)

        print(f"Shuffled file written: {output_filename}")

if __name__ == "__main__":
    main()