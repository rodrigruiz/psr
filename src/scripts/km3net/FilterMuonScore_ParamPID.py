import pandas as pd
import argparse
import numpy as np
from pathlib import Path

def filter_events(input_file, threshold, output_dir):
    """
    Load HDF5 file, filter events based on muon score, and save filtered data.
    """
    try:
        df = pd.read_hdf(input_file, "summary", dtype=np.float32)
    except KeyError:
        print(f"Error: {input_file} does not contain 'summary' table!")
        exit(1)

    # Check if 'muon_score' column exists
    if 'muon_score' not in df.columns:
        print(f"Warning: 'muon_score' column not found in {input_file}. No filtering applied.")
        df_filtered = df
    else:
        df_filtered = df[df["muon_score"] < threshold]

    # Generate filtered file name
    input_filename = Path(input_file).stem  # Extract filename without extension
    output_file = Path(output_dir) / f"filtered_{input_filename}.h5"

    # Save filtered data
    df_filtered.to_hdf(output_file, 'summary', format="table")
    print(f"Filtered data saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter HDF5 files based on muon score")
    parser.add_argument("--input", required=True, help="Path to input HDF5 file")
    parser.add_argument("--threshold", type=float, required=True, help="Muon score threshold")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory for filtered file")

    args = parser.parse_args()
    filter_events(args.input, args.threshold, args.output_dir)
