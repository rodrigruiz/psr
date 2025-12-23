"""
Add weights to the mc_table in HDF5 KM3NeT files based on a consolidated external weight file.

Usage: 
    AddWeightsKM3NeT_new.py -i INPUT_FILES... -w=<weights_file> -o OUTPUT_DIR [--detector=<detector>]

Options:
    -h --help                            Show this help message.
    -i --input_files INPUT_FILES         Input files or file pattern.
    -w --weights_file=<weights_file>     HDF5 file containing all weights.
    -o --output_dir OUTPUT_DIR           Output directory.
       --detector=<detector>             Detector name: 'arca' or 'orca' [default: arca]
"""

from docopt import docopt
import os, glob, re
import h5py
import numpy as np
import pandas as pd
from astropy.table import Column, Table
from astropy.io.misc.hdf5 import write_table_hdf5
from km3astro.io import load_hdf5_tables


def extract_run_id(filename: str) -> str:
    """Extract the last 8-digit run_id from the filename."""
    matches = re.findall(r'(?<!\d)\d{8}(?!\d)', filename)
    return matches[-1] if matches else None


def fetch_weights(mc_table, run_id, weight_table, detector):
    """Map weights to the mc_table using (run_id, event_id) keys."""
    # Filter for current file's run_id
    run_id = int(run_id)

    # Filter weight_table
    subset = weight_table[
        (weight_table['detector'] == detector) &
        (weight_table['run_id'] == run_id)
    ]

    # Check columns
    for col in ['event_id', 'new_weight']:
        if col not in subset.columns:
            raise ValueError(f"Missing required column: {col}")

    # Build event_id → new_weight map
    weight_map = {int(eid): float(w) for eid, w in zip(subset['event_id'], subset['new_weight'])}

    # Apply to mc_table
    weights = []
    run_ids_debug = []
    event_ids_debug = []

    for eid in mc_table['event_id']:
        weight = weight_map.get(int(eid), np.nan)
        weights.append(weight)

        run_ids_debug.append(int(run_id))
        event_ids_debug.append(int(eid))

    return weights, run_ids_debug, event_ids_debug


def write_all_tables(tables, output_file_path, weight_column, run_id_column, event_id_column):
    with h5py.File(output_file_path, 'w') as h5file:
        for attr_name in tables.__dict__:
            table = getattr(tables, attr_name)
            if table is None:
                continue

            if attr_name == 'mc_table':
                table.add_column(weight_column)
                table.add_column(run_id_column)
                table.add_column(event_id_column)
                group = h5file.require_group("MC")
                write_table_hdf5(table, group, path='MC_EVENTS', serialize_meta=True)
            elif attr_name == 'header_table':
                write_table_hdf5(table, h5file, path='HEADER', serialize_meta=True)
            elif attr_name == 'id_table':
                write_table_hdf5(table, h5file, path='ID', serialize_meta=True)
            elif attr_name in ['reco_table', 'fitinf_table']:
                group = h5file.require_group("RECO")
                path = "RECO_EVENTS" if attr_name == 'reco_table' else "FITINF"
                write_table_hdf5(table, group, path=path, serialize_meta=True)
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

    weights_file = str(data['weights_file'])
    weight_table = Table.read(weights_file, path='weights')

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    for file in input_files:
        folder_path, file_name = os.path.split(file)
        file_base = os.path.splitext(file_name)[0]
        detector = str(data['detector']).lower()

        # Remove _<detector>_classified if present
        suffix_to_remove = f"_{detector}_classified"
        if file_base.endswith(suffix_to_remove):
            file_base = file_base[: -len(suffix_to_remove)]

        # Extract run_id
        run_id = extract_run_id(file_name)
        if run_id is None:
            print(f"Could not extract run_id from {file_name}, skipping.")
            continue

        output_file_name = f"{file_base}_with_weights.h5"
        output_file_path = os.path.join(data['output_dir'], output_file_name)

        print(f"\nProcessing: {file}")
        print(f"Using run_id: {run_id}")
        tables = load_hdf5_tables(file)

        if not hasattr(tables, 'mc_table') or tables.mc_table is None:
            print(f"Warning: No mc_table found in {file}, skipping.")
            continue

        weights, run_ids_debug, event_ids_debug = fetch_weights(tables.mc_table, run_id, weight_table, detector)
        weight_column = Column(weights, name='atm_weight')
        run_id_column = Column(run_ids_debug, name='run_id_test')
        event_id_column = Column(event_ids_debug, name='event_id_test')

        write_all_tables(tables, output_file_path, weight_column, run_id_column, event_id_column)
        print(f"Written output file with weights: {output_file_path}")


if __name__ == "__main__":
    main()
