""" Combines multiple eventlist objects generated from KM3NeT data to a single one.

Usage: CombineEventListsKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR 

Options:
  -h --help                              Show this help message
  -i --input_files INPUT_FILES...        Input files
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
    arguments = docopt(__doc__)

    input_files = arguments['--input_files']
    output_dir = arguments['--output_dir']

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    CombinedEventList = None

    for file in input_files:
        with h5py.File(file, 'r') as h5_file:
            EventList = EL.readEventList(h5_file)

            if CombinedEventList is None:
                CombinedEventList = EventList
            else:
                CombinedEventList = vstack([CombinedEventList,EventList])
            
    CombinedEventList.sort("time")
    print("CombinedEventlist:")
    print(CombinedEventList)

    event_count = len(CombinedEventList)


    # Extract common prefix from input filenames
    common_prefix = os.path.commonprefix(input_files)
    # Remove any trailing non-alphanumeric characters from common prefix
    common_prefix = os.path.basename(common_prefix).rstrip("_-.")
    if not common_prefix:
        common_prefix = "Test"
    #output_file = os.path.join(output_dir, f"{common_prefix}_combined_eventlist.hdf5")
    output_file = os.path.join(output_dir, f"{common_prefix}_combined_eventlist_{event_count}events.hdf5")


    # Save the combined EventList
    CombinedEventList.write(output_file, format='hdf5', path='data', overwrite=True, serialize_meta = True)
    print(f"Combined EventList saved to {output_file}")

    
if __name__ == "__main__":
    main()
    