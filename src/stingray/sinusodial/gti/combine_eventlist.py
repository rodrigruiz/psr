#!/usr/bin/env python3
"""
Combine timeseries to one timeseries.

Usage:
  combine_eventlist.py --directory=<directory>
  combine_eventlist.py (-h | --help)

Options:
  -h --help                   Show this help message and exit.
  --directory=<directory>     Directory where timeseries files are stored.
"""
import os, glob
from docopt import docopt
from astropy.timeseries import TimeSeries
from astropy.table import vstack
#from analyse_timeseries import combine_timeseries

def combine_eventlist(directory):
    file_pattern = os.path.join(directory, 'eventlist_*.dat')
    timeseries_files = glob.glob(file_pattern)

    if not timeseries_files:
        print(f"No timeseries files found in the directory: {directory}")
        return

    # Sort files by the integer value i in the filename
    timeseries_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
    
    combined_ts = None
    
    for file_path in timeseries_files:
        ts = TimeSeries.read(file_path, format='ascii', time_column='time', time_format='mjd')
        #print( combined_ts )
        if combined_ts is None:
            combined_ts = ts
        else:
            combined_ts = vstack([combined_ts, ts])
    combined_ts.write('combined_eventlist.dat', format='ascii.ecsv', overwrite=True) 
    
if __name__ == '__main__':
    arguments = docopt(__doc__)
    combine_eventlist(
        arguments['--directory']
    )