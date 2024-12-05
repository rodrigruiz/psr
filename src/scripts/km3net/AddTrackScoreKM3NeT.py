"""
Extract the track_score column from h5 km3net files that were classified using ParamPID.
The filename of the rootfilename should be given, as well as the folder where the classified files are located.
Those are named rootfilename_scored.h5 and this script expects this to be the case because it appends the _scored.h5 part to the filepath of the rootfile.

Usage: AddTrackScoreKM3NeT.py -f input_folder -i INPUT_FILES... -o OUTPUT_DIR [--detector=<detector>]

Options:
  -h --help                              Help
  -f --input_folder INPUT_FOLDER         Folder of the input files
  -i --input_files INPUT_FILES           Input filenames
  -o --output_dir OUTPUT_DIR             Output directory
     --detector=<detector>               Name of the Detector. 'arca' or 'orca' [default: arca]
"""

from docopt import docopt
import os
import glob
import h5py
from km3astro import io


def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]

    input_folder = data['input_folder']
    input_files = []
    for pattern in data['input_files']:
        full_pattern = os.path.join(input_folder, pattern)
        input_files.extend(glob.glob(full_pattern))
    input_files.sort()

    det_name = str(data['detector'])

    if not input_files:
        print(f"No files matching pattern in folder {input_folder}: {data['input_files']}")
        return

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    for file in input_files:
        if not os.path.exists(file):
            print(f"Warning: File {file} does not exist. Skipping.")
            continue

        folder_path, file_name = os.path.split(file)
        file_base, _ = os.path.splitext(file_name)
        output_filename = os.path.join(data['output_dir'], f"{file_base}_scored.h5")
        
        print(f"Processing: {file}")
        print(f"Output: {output_filename}")

        # Convert root file to HDF5
        #io.root_to_hdf5(file, output_file=output_filename)


if __name__ == "__main__":
    main()