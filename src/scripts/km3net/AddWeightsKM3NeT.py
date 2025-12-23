"""
Add weights to the mc_table in HDF5 KM3NeT files based on external weight files.
The weight files are expected to be named like the input file with '_atm_neutrino_weights.hdf5' appended before the extension.

Usage: 
    AddWeightsKM3NeT.py -i INPUT_FILES... -p=<weights_folder> -o OUTPUT_DIR [--detector=<detector>]

Options:
    -h --help                                      Show this help message.
    -i --input_files INPUT_FILES                   Input files or file pattern.
    -p --weights_folder=<weights_folder>           Folder containing the weight files.
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


def calculate_weights(mc_table, weight_filepath):
    """Load weight values and map them to the mc_table based on event_id."""
    try:
        df = pd.read_hdf(weight_filepath)
    except Exception as e:
        raise IOError(f"Error reading weight file: {weight_filepath}. Details: {e}")

    if 'event_id' not in df.columns or 'new_weight' not in df.columns:
        raise ValueError(f"Weight file {weight_filepath} is missing required columns 'event_id', 'new_weight'.")

    weight_map = dict(zip(df['event_id'], df['new_weight']))
    event_ids = mc_table['event_id']

    unmatched = [eid for eid in event_ids if eid not in weight_map]
    if unmatched:
        print(f"Warning: {len(unmatched)} unmatched event_ids. First 10: {unmatched[:10]}")

    return [weight_map.get(eid, np.nan) for eid in event_ids]


def write_all_tables(tables, output_file_path, weight_column):
    with h5py.File(output_file_path, 'w') as h5file:
        for attr_name in tables.__dict__:
            table = getattr(tables, attr_name)
            if table is None:
                continue

            if attr_name == 'header_table':
                write_table_hdf5(table, h5file, path='HEADER', serialize_meta=True)
            elif attr_name == 'id_table':
                write_table_hdf5(table, h5file, path='ID', serialize_meta=True)
            elif attr_name in ['reco_table', 'fitinf_table']:
                group = h5file.require_group("RECO")
                path = "RECO_EVENTS" if attr_name == 'reco_table' else "FITINF"
                write_table_hdf5(table, group, path=path, serialize_meta=True)
            elif attr_name == 'mc_table':
                table.add_column(weight_column)
                group = h5file.require_group("MC")
                write_table_hdf5(table, group, path='MC_EVENTS', serialize_meta=True)
            else:
                write_table_hdf5(table, h5file, path=attr_name.upper(), serialize_meta=True)


def main():
    arguments = docopt(__doc__)
    data = {key.replace("-", ""): arguments[key] for key in arguments}

    input_files = []
    for pattern in data['input_files']:
        input_files.extend(glob.glob(pattern))
    input_files.sort()

    if not input_files:
        print(f"No files matching pattern: {data['input_files']}")
        return

    weights_folder = str(data['weights_folder'])

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    for file in input_files:
        folder_path, file_name = os.path.split(file)
        file_base = os.path.splitext(file_name)[0]

        # Get detector name as string
        detector = str(data['detector']).lower()

        # Remove the '_<detector>_classified' suffix if it exists
        suffix_to_remove = f"_{detector}_classified"
        if file_base.endswith(suffix_to_remove):
            file_base = file_base[: -len(suffix_to_remove)]

        weight_file_name = f"{file_base}_atm_neutrino_weights.hdf5"
        weight_filepath = os.path.join(weights_folder, weight_file_name)

        output_file_name = f"{file_base}_with_weights.h5"
        output_file_path = os.path.join(data['output_dir'], output_file_name)

        print(f"Processing: {file}")
        tables = load_hdf5_tables(file)
        print(f"Available tables: {list(tables.__dict__.keys())}")

        if not hasattr(tables, 'mc_table') or tables.mc_table is None:
            print(f"Warning: No mc_table found in {file}, skipping.")
            continue

        print(f"Reading weights from: {weight_filepath}")
        weights = calculate_weights(tables.mc_table, weight_filepath)
        weight_column = Column(weights, name='atm_weight')

        write_all_tables(tables, output_file_path, weight_column)
        print(f"Written output file with weights: {output_file_path}")


if __name__ == "__main__":
    main()
