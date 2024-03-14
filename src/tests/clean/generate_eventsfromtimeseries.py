#!/usr/bin/env python3
"""
Generate event times from a timeseries.

Usage:
  generate_eventsfromtimeseries.py [--directory=<directory> | --timeseries=<timeseries>] [--output=<outputfilepath>] [--data_column_name=<data_column_name>] [--time_column_name=<time_column_name> --data_column_name=<data_column_name> --format=<format>]
  generate_eventsfromtimeseries.py (-h | --help)

Options:
  -h --help                   Show this help message and exit.
  --directory=<directory>     Directory where timeseries files are stored.
  --time_column_name=<time_column_name> Time column name [default:time]
  --data_column_name=<data_column_name> data column name [default:counts]
"""

import os, glob
from collections.abc import Iterable
import numpy as np
from astropy.timeseries import TimeSeries, BinnedTimeSeries
from astropy.table import Table
from astropy.time import Time
from docopt import docopt

def load_timeseries(directory, input_file='timeseries_*.dat'):
    """Loads astropy.timeseries.TimeSeries from a single or multiple files.
    
    Parameters
    ----------
        directory : str
            Directory containing the 'timeseries_i.dat' files (i is an integer value starting from 0).
            
    Returns
    -------
        list or astropy.timeseries.TimeSeries object
            List of astropy.timeseries.TimeSeries objects or a single astropy.timeseries.TimeSeries object.
        
    """
    
    #file_pattern = os.path.join(directory, 'timeseries_*.dat')
    file_pattern = os.path.join(directory, input_file)
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
    
def create_eventlist(ts=None, directory=None, output='eventlist', time_column_name='time', data_column_name='counts', scaling=1, format='ascii'):
    """Creates an eventlist from an astropy.timeseries.TimeSeries.
    Assumes evenly distributed events inside a timebin.
    
    Parameters
    ----------
        ts : astropy.timeseries.TimeSeries
            TimeSeries to create eventlist from.
            `ts` and `direcotry` are mutually exclusive.
            
        directory : str
            Directory containing the 'timeseries_i.dat' files (i is an integer value starting from 0).
            `ts` and `direcotry` are mutually exclusive.
    
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
        Saves the eventlist to a file in the ascii.ecsv or hdf5 format.

    """
    np.random.seed()
    
    if ts is None:
        ts = load_timeseries(directory)
        print(ts)
    elif directory is None:
        ts = ts
    else:
        raise ValueError(
            "ts or directory is not specified.")
    
    print(ts)
    if isinstance(ts, Iterable):
        for i, s in enumerate(ts):
            output_file = output + f'_{i}'
            _generate_events(s, output=output_file, time_column_name=time_column_name, data_column_name=data_column_name, format=format)
    else:
        _generate_events(ts, output=output, time_column_name=time_column_name, data_column_name=data_column_name, scaling=scaling, format=format)
    
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
    
def _generate_events(ts, output='eventlist', time_column_name='time', data_column_name='counts', scaling=1, format='ascii'):
    """Helping function.
       Generates events and writes the eventlist to a ascii.escv file.
    """
    
    event_list = None
    for i, bins in enumerate(ts[:-1]):
        if isinstance(ts, TimeSeries):
            events = np.random.uniform( bins[time_column_name].value, ts[i+1][time_column_name].value, int(bins[data_column_name]/scaling) )
        elif isinstance(ts, BinnedTimeSeries):
            events = np.random.uniform( bins[time_column_name].value, ts[time_column_name][i+1].value, int(bins[data_column_name].value/scaling) )
        else:
            raise ValueError(
                "ts should be an astropy.timeseries.TimeSeries or an astropy.timeseries.BinnedTimeSeries."
                            )
            
        events = np.sort(events)
        #print(events)
        if event_list is None:
            event_list = events
        else:
            event_list = np.concatenate((event_list, events))

    #event_list = Table([event_list], names=[time_column_name])
    print(event_list)
    event_list = Table([event_list], names=['time'])
    
    if format == 'hdf5':
        event_list.write(output + '.hdf5', format='hdf5', overwrite=True)
    else:
        event_list.write(output + '.dat', format='ascii.ecsv', overwrite=True)
    
if __name__ == '__main__':
    arguments = docopt(__doc__)
    directory = arguments['--directory']
    timeseries = arguments['--timeseries']
    output = arguments['--output']
    
    if directory is None:
        if output is None:
            create_eventlist(ts=timeseries, time_column_name=arguments['--time_column_name'], data_column_name=arguments['--data_column_name'], format=arguments['--format'])
        else:
            create_eventlist(ts=timeseries, output=output, time_column_name=arguments['--time_column_name'], data_column_name=arguments['--data_column_name'], format=arguments['--format'])
            
    else:
        if output is None:
            create_eventlist(directory=directory)
        else:
            create_eventlist(direcotry=directory, output=output)
    
    
    