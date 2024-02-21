#!/usr/bin/env python3
"""
Generate event times from a timeseries.

Usage:
  generate_eventsfromtimeseries.py --directory=<directory>
  generate_eventsfromtimeseries.py (-h | --help)

Options:
  -h --help                   Show this help message and exit.
  --directory=<directory>     Directory where timeseries files are stored.
"""

import numpy as np
from astropy.timeseries import TimeSeries
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
from docopt import docopt
import matplotlib.pyplot as plt
from analyse_timeseries import load_timeseries
from astropy.io import fits

if __name__ == '__main__':
    arguments = docopt(__doc__)
    directory = arguments['--directory']
    
    ts = load_timeseries(directory)
    #print(ts[0])
    np.random.seed()
    
    event_list = None
    for i, bins in enumerate(ts[:-1]):
        #print(bins['time'].value,  ts[i+1]['time'].value, bins['flux'])
        events = np.random.uniform( bins['time'].value, ts[i+1]['time'].value, int(bins['countrate']) )
        np.sort(events)
        #print(len(events))
        if event_list is None:
            event_list = events
        else:
            event_list = np.concatenate((event_list, events))
    
    #print(event_list, np.shape(event_list))
    #event_list = TimeSeries(time=event_list, format='unix')
    event_list = Table([event_list], names=['time'])
    file_name = 'eventlist.dat'
    event_list.write(file_name, format='ascii.ecsv', overwrite=True)
    #fits.writeto(file_name, event_list)
    #np.savetxt('eventlist.txt', event_list)
    #with open("eventlist.txt", "w") as o:
    #    for event in event_list:
    #        print("{}".format(event), file=o) 
    