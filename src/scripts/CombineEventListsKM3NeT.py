""" Combines multiple eventlist objects generated from KM3NeT data to a single one.

Usage: CombineEventListsKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR 

Options:
  -h --help                              Show this help message
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output directory  
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
    # Getting Key-Argument-Pairs that are passed to the script
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

    # Construct TimeSeries for each run
    CombinedEventList = None

    for file in input_files:
        with h5py.File(file) as h5_file:
            EventList = EL.readEventList(h5_file)

            if CombinedEventList is None:
                CombinedEventList = EventList
            else:
                CombinedEventList = vstack([CombinedEventList,EventList])
            
    CombinedEventList.sort("time")
    print("CombinedEventlist:")
    print(CombinedEventList)

    # Extract common prefix from input filenames
    common_prefix = os.path.commonprefix(input_files)
    # Remove any trailing non-alphanumeric characters from common prefix
    common_prefix = os.path.basename(common_prefix).rstrip("_-.")
    if not common_prefix:
        common_prefix = "Test"
    output_file = os.path.join(data['output_dir'], f"{common_prefix}_combined_eventlist.h5")


    # Save the combined EventList
    CombinedEventList.write(output_file, format='hdf5', path='data', overwrite=True, serialize_meta = True)
    print(f"Combined EventList saved to {output_file}")

    
if __name__ == "__main__":
    main()
    