""" Split Time Series

Usage: split_timeseries.py -i INPUT_FILES -o OUTPUT_DIR

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output directory

"""
import h5py
import numpy as np
import matplotlib.pyplot as plt
import plens.TimeSeries as TS
#import scripts 
import sys
sys.path.insert(0, '~software/Psr/src/stingray/clean')
from plot_timeseries import plot_timeseries
from generate_split_timeseries import split_timeseries
from generate_eventsfromtimeseries import create_eventlist

plt.style.use('~/software/Psr/src/stingray/clean/latex.mplstyle')

# do analysis on 0-ending runs only!!!
# Read TimeSeries
#if __name__ == "__main__":
    
    #file = '/home/wecapstor3/capn/mppi19/ANTARES/data/rates/hdf5/Antares_051870_total_rates.hdf5'
    
    
file = 'Antares_051870_total_rates_combined.hdf5'
with h5py.File(file) as input_file:
    ts, timeslice_duration = TS.readTimeSeries(input_file)
    split_timeseries(ts, len(ts), len(ts)/100, './data/051870/Antares_051870_total_rates', format='hdf5', antares=True, timeslice_duration=timeslice_duration)

"""
    ts['time_bin_start'].format = 'mjd'
    print(ts[:2])
    #print(TimeSeries)
    fig, ax = plt.subplots()
    plt.figure(figsize=(10, 5))
    ax.set_xticks(ts['time_bin_start'].value[::30000])
    plt.plot(ts['time_bin_start'].value[1:90000], ts['rateOn'][1:90000], color='blue')
    plt.xlabel('Time [MJD]')
    plt.ylabel('Counts [a.u.]')
    #plt.title(title)
    plt.grid(True)
    plt.savefig( 'Antares_ts.png', bbox_inches='tight' )
    #plt.legend()
    #plot_timeseries(file, output='Antares_ts', time_column_name='time_bin_start', data_column_name='rateOff', color='blue', label='countrates [Hz]', title='Antares run 051870 -- total rates Offshore')

'''
with h5py.File(file, 'r') as f:
    print("Keys: %s" % f.keys())
    #print(f['lcmID_dict'])
    print(f['header'].keys())
    print(f['total_rates'].attrs.keys())
'''
    