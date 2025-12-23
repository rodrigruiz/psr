"""  Check for gaps in the data. Determine good time intervals

Usage: FindGTIsKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR [--df=<df>] [--combine]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output directory 
     --df=<float>                        Time Resolution of the binned light curve [default: 60]
     --combine                           Whether to combine all the gtis or make separate gti files per input file
"""
#python3 psr/src/scripts/CorrectEventListKM3NeT.py -i '/home/hpc/capn/capn107h/software/eventlistTestOutput/mcv*' -o correctedeventlistTestOutput/ -s hdf5SourceFiles/Vela_X-1.h5

from docopt import docopt
import os, glob
import h5py as h5py
import re
import pickle
import matplotlib.pyplot as plt
import numpy as np
import plens.TimeSeries as TS
import plens.EventList as EL
from astropy.table import Table
from epochfolding.gtis import findGTIs, saveGTIs
from stingray import Lightcurve
import astropy.units as u



def main():
    # Getting Key-Argument-Pairs that are passed to the script
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    input_files = arguments['--input_files']
    df = float(arguments['--df'])

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])
    
    if data['combine']:
        gtis = []
        absolute_gtis = None

    for file in input_files:
        # Fetching filename for usage in output filename 
        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]

        output_file = data['output_dir'] + file_name + '_gtis'
        print(output_file, os.path.exists(output_file))

        # Read TimeSeries file and apply corrections using timing and source information
        with h5py.File(file) as input_file:

            EventList = EL.readEventList(input_file)
            print(EventList['time'])
            print(EventList['time'].value)
            tstart = EventList['time'][0].value
            print(tstart)
            print(df)
            lc = Lightcurve.make_lightcurve(EventList['time'], dt=df)
            print(lc)
            
            
            if data['combine']:
                gtis_file = findGTIs(lc.counts, lc.time) 
                absolute_gtis_file = gtis_file + tstart
                gtis += gtis_file
                print(absolute_gtis)
                print(absolute_gtis_file)
                print(len(absolute_gtis_file))
                if absolute_gtis is None:
                    absolute_gtis = absolute_gtis_file
                else:
                    absolute_gtis = np.concatenate([absolute_gtis,absolute_gtis_file])

            else:
                gtis = findGTIs(lc.counts, lc.time)
                print(gtis)
                absolute_gtis = gtis + tstart
                print(absolute_gtis)
                saveGTIs(absolute_gtis, output_file)



        if os.path.exists(output_file):
            os.remove(output_file)  # Remove the file if it already exists

    if data['combine']:
        combined_output_file = data['output_dir'] + 'gtifile_combined_gtis'
        absolute_gtis.sort()

        plt.figure(figsize=(10, 6))
        plt.plot(lc.time, lc.counts, label="Counts", color="blue")
        for start, end in gtis:
            plt.axvline(start, linestyle="dotted", color="gray", alpha=0.8)
            plt.axvline(end, linestyle="dotted", color="gray", alpha=0.8)
            plt.axvspan(start, end, color="gray", alpha=0.2) 
        plt.xlabel("Time (s)", fontsize=14)
        plt.ylabel("Counts", fontsize=14)
        plt.title("Lightcurve with highlighted GTIs", fontsize=16)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("Testplot_lightcurve_withGTIs.png")
        
        print(len(absolute_gtis))
        print(absolute_gtis)
        saveGTIs(absolute_gtis, combined_output_file)

    
if __name__ == "__main__":
    main()
    