"""
Usage: TestScript.py INPUT_FILES... -o OUTPUT_DIR

Options:
  INPUT_FILES                      Input files
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

    print(data)
    input_files = []
    for pattern in data['INPUT_FILES']:
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
    