"""
Usage:
  TestScript.py -i INPUT_FILES... -o OUTPUT_DIR

Options:
  -i --input_files INPUT_FILES    Input file paths.
  -o --output_dir OUTPUT_DIR      Output directory.
"""

# python3 psr/src/scripts/CombineEventListsKM3NeT.py -i "/home/hpc/capn/capn107h/software/correctedeventlistTestOutput/*" -o combinedeventlistTestOutput/


from docopt import docopt
import os
import glob
import h5py
import numpy as np
from astropy.table import vstack, Table
import plens.EventList as EL

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

    output_dir = data['output_dir']

    print(output_dir)
    
    for file in input_files:
        print(file)

    
if __name__ == "__main__":
    main()
    