#!/usr/bin/env python3
"""
Generate event times from a timeseries.

Usage:
  generate_eventsfromtimeseries.py --directory=<directory> [--output=<outputfilepath>] [--data_column_name=<data_column_name>]
  generate_eventsfromtimeseries.py (-h | --help)

Options:
  -h --help                   Show this help message and exit.
  --directory=<directory>     Directory where timeseries files are stored.
"""

import os, glob
from collections.abc import Iterable
import numpy as np
from astropy.timeseries import TimeSeries
from astropy.table import Table
from astropy.time import Time
from docopt import docopt

def load_timeseries(directory):
    """Loads astropy.timeseries.TimeSeries from multiple files.
    
    Parameters
    ----------
        directory : str
            Directory containing the 'timeseries_i.dat' files (i is an integer value starting from 0).
            
    Returns
    -------
        list or astropy.timeseries.TimeSeries object
            List of astropy.timeseries.TimeSeries objects or a single astropy.timeseries.TimeSeries object.
        
    """
    
    file_pattern = os.path.join(directory, 'timeseries_*.dat')
    timeseries_files = glob.glob(file_pattern)

    if not timeseries_files:
        print(f"No timeseries files found in the directory: {directory}")
        return

    # Sort files by the integer value i in the filename
    timeseries_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
    
    if len(timeseries_files) == 1:
        return TimeSeries.read(timeseries_files[0], format='ascii', time_column='time', time_format='mjd')
    else:
        timeseries = []
        for file in timeseries_files:
            ts = TimeSeries.read(file, format='ascii', time_column='time', time_format='mjd')
            timeseries.append(ts)
        return timeseries
    
def create_eventlist(directory, output='eventlist', data_column_name='counts'):
    """Creates an eventlist from an astropy.timeseries.TimeSeries.
    Assumes evenly distributed events inside a timebin.
    
    Parameters
    ----------
        directory : str
            Directory containing the 'timeseries_i.dat' files (i is an integer value starting from 0).
    
    Optional Parameters
    -------------------
        output : str
            Output filename or entire path to the outputfile (without the fileending!).
            Defaults to 'timeseries'.
        
        data_column_name : str
            Column name of the data to plot.
            Defaults to 'counts'.
            
    Returns
    -------
        Saves the eventlist to a file in the ascii.ecsv format.

    """
    np.random.seed()
    
    ts = load_timeseries(directory)
    
    if isinstance(ts, Iterable):
        for i, s in enumerate(ts):
            output_file = output + f'_{i}'
            _generate_events(s, output=output_file)
    else:
        _generate_events(ts)
    
    '''
    event_list = None
    for i, bins in enumerate(ts[:-1]):
        events = np.random.uniform( bins['time'].value, ts[i+1]['time'].value, int(bins[data_column_name]) )
        events = np.sort(events)
        
        if event_list is None:
            event_list = events
        else:
            event_list = np.concatenate((event_list, events))

    event_list = Table([event_list], names=['time'])
    event_list.write(output + '.dat', format='ascii.ecsv', overwrite=True)
    '''
    
def _generate_events(ts, output='eventlist', data_column_name='counts'):
    """Helping function.
       Generates events and writes the eventlist to a ascii.escv file.
    """
    
    event_list = None
    for i, bins in enumerate(ts[:-1]):
        events = np.random.uniform( bins['time'].value, ts[i+1]['time'].value, int(bins[data_column_name]) )
        events = np.sort(events)
        
        if event_list is None:
            event_list = events
        else:
            event_list = np.concatenate((event_list, events))

    event_list = Table([event_list], names=['time'])
    event_list.write(output + '.dat', format='ascii.ecsv', overwrite=True)
    
if __name__ == '__main__':
    arguments = docopt(__doc__)
    directory = arguments['--directory']
    output = arguments['--output']
    
    if output is None:
        create_eventlist(directory)
    else:
        create_eventlist(directory, output=output)
    
    
    