""" Load KM3NeT root data files and convert them to astropy tables. 

Usage: ConvertFiles.py -i INPUT_FILES... -o OUTPUT_DIR

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files or file pattern
  -o --output_dir OUTPUT_DIR             Output directory

"""
#python3 ConvertFiles.py -i '/home/hpc/capn/mppi104h/wecapstor3/out/ARCA/KM3NeT_00000133/v8.1/reco/*.root' -o TestOutput/

from docopt import docopt
import os, glob
import re
import h5py as h5py
from km3astro import io


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
        output_filename = data['output_dir'] + file_name + ".h5"
        io.root_to_hdf5(file, output_file = output_filename)

            
if __name__ == "__main__":
    main()