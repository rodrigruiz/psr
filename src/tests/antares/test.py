import h5py
import numpy as np
from docopt import docopt
import os, glob
import re
import plens.TimeSeries as TS
from astropy.time import Time, TimeDelta

"""
Porb = 2.0869953 #day
T_pi2 = 2455073.68504 - 2400000.5 # jd
Tnod = T_pi2 + Porb/2 #mjd
t = Time(Tnod, format='mjd').tdb
print(Tnod, t, t.value*86400, t.unix, type(t.unix))
# test 1
file = './data/Antares_051870_total_rates_0.hdf5'
with h5py.File(file) as h5_file:
    ts, timeslice_duration = TS.readTimeSeries(h5_file)
    print(ts['time_bin_start'][0], ts['time_bin_start'].format, ts['time_bin_start'][0].mjd, ts['time_bin_start'][0].tdb, ts['time_bin_start'][0].mjd * 86400)

file = './data/Antares_051870_total_rates_0_corrected.hdf5'
with h5py.File(file) as h5_file:
    ts, timeslice_duration = TS.readTimeSeries(h5_file)
    print(ts['time_bin_start'][0])
"""

# test 2 
files = glob.glob('./data/Antares_05*0_total_rates_combined.hdf5')
print(files)
for file in files:
    with h5py.File(file) as f:
        print(len(f['timeseries/time_bin_start']))
    #new = np.array_split(f['timeseries/time_bin_start'], np.ceil(len(f['timeseries/time_bin_start'])/100))
    #print(len(new))
    #print(new)
"""