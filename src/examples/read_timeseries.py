#!/usr/bin/env python3
"""
Analyze time series data and create a combined plot of folded time series.

Usage:
  analyze_timeseries.py --directory=<directory> --period=<period> --reference_mjd=<reference_mjd>
  analyze_timeseries.py (-h | --help)

Options:
  -h --help                   Show this help message and exit.
  --directory=<directory>     Directory where timeseries files are stored.
  --period=<period>           Period for folding.
  --reference_mjd=<reference_mjd> Reference MJD for folding.
"""

import os
import glob
import numpy as np
from astropy.timeseries import TimeSeries
import matplotlib.pyplot as plt
from docopt import docopt
from astropy import units as u  # Added for Quantity
from astropy.time import Time  # Added for Time
from astropy.timeseries import aggregate_downsample
from astropy.table import vstack

# chi2 test to check against a constant profile
def get_chi2(unfolded_ts, folded_ts, period, flux_column_name='flux'):
    # std of unfolded time series - load the data from the different files first, combine and then calculate std
    # x_bar: from Larsson(1996) mean of unfolded data since we want to test against a constant profile (?folded time series)
    # std: from Larsson(1996) std of unfolded data divided by #points in bin = test period / bin width
    
    # mean and std of unfolded data
    mean, std = get_stats_timeseries( unfolded_ts )
    chi2 = 0
    # total number of bins is given by the length of the folded data
    print(f'length of folded timeseries:{len(folded_ts)}')
    for i in folded_ts[flux_column_name]:
        # number of data points in bin i - for this case is the same in every bin, but could be different in experimental data if the resolution changes or the are gaps in the data
        #print(f'number of data points in bin: {n_i}')
        n_i = ((unfolded_ts['time'].value[-1]) - (unfolded_ts['time'].value[0])) / period
        print(f'number of data points in bin: {n_i}')
        chi2 += ( i - mean )**2 / ( std / np.sqrt(n_i) )**2
    return chi2

def calculate_stats(ts):
     # Assuming you have a TimeSeries object named 'ts' with a 'flux' column

     # Calculate the length of the 'flux' column
     num_rows = len(ts)

     # Calculate the split points
     first_third_start = 0
     first_third_end = num_rows // 3
     second_third_start = first_third_end
     second_third_end = 2 * (num_rows // 3)
     last_third_start = second_third_end
     last_third_end = num_rows

     # Divide the data into three equal parts
     first_third = ts[first_third_start:first_third_end]
     second_third = ts[second_third_start:second_third_end]
     last_third = ts[last_third_start:last_third_end]

     # Calculate the mean and standard deviation for the first third
     first_third_mean = np.mean(first_third['flux'])
     first_third_std = np.std(first_third['flux'])

     # Calculate the mean and standard deviation for the second third
     second_third_mean = np.mean(second_third['flux'])
     second_third_std = np.std(second_third['flux'])

     # Calculate the mean and standard deviation for the last third
     last_third_mean = np.mean(last_third['flux'])
     last_third_std = np.std(last_third['flux'])

     print("First Third Mean:", first_third_mean)
     print("First Third Standard Deviation:", first_third_std)
     print("Second Third Mean:", second_third_mean)
     print("Second Third Standard Deviation:", second_third_std)
     print("Last Third Mean:", last_third_mean)
     print("Last Third Standard Deviation:", last_third_std)

def calculate_segment_stats(time_series, segment_length=10, flux_column_name='flux'):
    # Extract the 'flux' column from the TimeSeries
    flux = time_series[flux_column_name]

    # Convert the TimeSeries times to Julian Date (jd)
    times_jd = np.sort(time_series.time.jd)
   # Calculate the total number of data points
    num_data_points = len(time_series)

    # Calculate the number of segments
    num_segments = num_data_points // segment_length

    # Initialize arrays to store the results
    segment_times = []
    mean_flux_values = []
    std_flux_values = []

    print(f"the number of segments is {num_segments}")

    for i in range(num_segments):
        # Calculate the start and end indices for the current segment
        start_idx = i * segment_length
        end_idx = (i + 1) * segment_length

        if i==int(num_segments / 2):
            print(f"start_idx={start_idx}, end_idx={end_idx}")

        # Extract the data within the current segment
        segment_flux = flux[start_idx:end_idx]
        segment_times_jd = times_jd[start_idx:end_idx]
        # Calculate the central time of the segment (in Julian Date)
        # segment_time_jd = np.mean(segment_times_jd)
        segment_time_jd = times_jd[start_idx+int(segment_length/2)]
        if i==int(num_segments/2):
            print(f"segment_time = {segment_time_jd}")

        # Convert the central time back to a Time object if needed
        #segment_time = time_series.time[0].from_jd(segment_time_jd)

        # Calculate the mean and standard deviation of the segment's flux
        segment_mean_flux = np.mean(segment_flux)
        segment_std_flux = np.std(segment_flux)
        if i==int(num_segments/2):
            print(segment_flux)
            print(f'mean flux = {segment_mean_flux}')
        # Append the results to the respective arrays
        segment_times.append(segment_time_jd)
        mean_flux_values.append(segment_mean_flux)
        std_flux_values.append(segment_std_flux)

    return np.array(segment_times), np.array(mean_flux_values), np.array(std_flux_values)



# def calculate_segment_stats(time_series, n):
#     # Extract the 'times' and 'flux' columns from the TimeSeries
#     # times = time_series['times']
#     # Convert the TimeSeries times to Julian Date (jd)
#     times = time_series.time.jd

#     flux = time_series['flux']

#     # Calculate the number of data points in each segment
#     segment_length = len(time_series) // n

#     # Initialize arrays to store the results
#     segment_times = []
#     mean_flux_values = []
#     std_flux_values = []

#     for i in range(n):
#         # Calculate the start and end indices for the current segment
#         start_idx = i * segment_length
#         end_idx = (i + 1) * segment_length if i < n - 1 else len(time_series)

#         # Extract the data within the current segment
#         segment_flux = flux[start_idx:end_idx]

#         # Calculate the central time of the segment
#         segment_time = np.mean(times[start_idx:end_idx])

#         # Calculate the mean and standard deviation of the segment's flux
#         segment_mean_flux = np.mean(segment_flux)
#         segment_std_flux = np.std(segment_flux)

#         # Append the results to the respective arrays
#         segment_times.append(segment_time)
#         mean_flux_values.append(segment_mean_flux)
#         std_flux_values.append(segment_std_flux)

#     return np.array(segment_times), np.array(mean_flux_values), np.array(std_flux_values)

def combine_timeseries(directory, flux_column_name='flux'):
    '''combines the timeseries from different files into a single timeseries.'''
    
    file_pattern = os.path.join(directory, 'timeseries_*.dat')
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
        
    return combined_ts

def get_stats_timeseries(timeseries, column_name='flux'):
    '''calculates mean and std of a timeseries column.'''
    flux = timeseries[column_name]
    return np.mean(flux), np.std(flux)

def analyze_split_timeseries(directory, period, reference_mjd):
    # Find all files named 'timeseries_*.dat' in the specified directory
    file_pattern = os.path.join(directory, 'timeseries_*.dat')
    timeseries_files = glob.glob(file_pattern)

    if not timeseries_files:
        print(f"No timeseries files found in the directory: {directory}")
        return

    # Sort files by the integer value i in the filename
    timeseries_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))

    # Initialize an empty TimeSeries to store the combined folded data
    combined_folded_ts = None
    
    # Create a combined plot
    plt.figure(figsize=(10, 6))
    for file_path in timeseries_files:
        ts = TimeSeries.read(file_path, format='ascii', time_column='time', time_format='mjd')
        print(f'unfolded timeseries: {ts}')
        print(f'length of unfolded timeseries: {len(ts)}')
        # Plot the original time series
        #plt.plot(ts.time.jd, ts['flux'], '.', markersize=1, label=file_path)

        # Fold the time series with the specified period and reference MJD
        
        epoch_time = Time(reference_mjd, format='mjd')  # Create a Time object for the reference MJD
        folded_ts = ts.fold(period=period * u.day, epoch_time=epoch_time)
        #plt.plot(folded_ts.time.jd, folded_ts['flux'], '.', markersize=1, color='red')
        print(f'folded timeseries: {folded_ts}')
        print( f'length of folded timeseries: {len(folded_ts)}' )
        # Sum the folded time series to the combined_folded_ts
        if combined_folded_ts is None:
            combined_folded_ts = folded_ts
        else:
            combined_folded_ts['flux'] += folded_ts['flux']
    
    plt.plot(combined_folded_ts.time.jd, combined_folded_ts['flux'], '.', markersize=1, color='red')
    # calculate chi2
    combined_unfolded_ts = combine_timeseries(directory)
    #print(float(period))
    chi2 = get_chi2( combined_unfolded_ts, combined_folded_ts, period)
    '''
    # Plot the folded time series
    stats = calculate_segment_stats(combined_folded_ts, 10)
    #print(stats[0])
    folded_ts.time.format = 'jd'
    plt.xlabel('Folded Time (JD)')
    plt.ylabel('flux [a.u.]')
    plt.title("Original Time Series")
    plt.figure(figsize=(10,6))
    plt.plot(combined_folded_ts.time.jd, combined_folded_ts['flux'], '.', markersize=1, color='black', label=file_path)
    ts_binned = aggregate_downsample(combined_folded_ts, time_bin_size=0.001 * u.day)
    ts_binned_2 = aggregate_downsample(combined_folded_ts, time_bin_size=0.001 * u.day, aggregate_func=np.std)
    #print(f'mean folded signal: {ts_binned}')
    #print(f'std of folded signal: {ts_binned_2}')
    plt.plot(ts_binned.time_bin_start.jd, ts_binned['flux'], 'r-', drawstyle='steps-post')
    plt.plot(ts_binned_2.time_bin_start.jd, ts_binned['flux'] + ts_binned_2['flux'], '-', color='orange', drawstyle='steps-post')
    plt.plot(ts_binned_2.time_bin_start.jd, ts_binned['flux'] - ts_binned_2['flux'], '-', color='orange', drawstyle='steps-post')
    #plt.plot(stats[0], stats[1], linestyle='None',marker='.', color='red')
    '''
    plt.xlabel('Folded Time (JD)')
    plt.ylabel('flux [a.u.]')
    
    #plt.legend(loc='upper right')
    plt.title(f'Folded Time Series. Period = {period} days.')
    plt.savefig('folded_timeseries.pdf')
    plt.show()
    
    return chi2

def analyze_timeseries(directory, period, reference_mjd):
    # calculate chi2
    ts = combine_timeseries(directory)
    
    # Create a combined plot
    plt.figure(figsize=(10, 6))
    
    print(f'unfolded timeseries: {ts}')
    print(f'length of unfolded timeseries: {len(ts)}')
    # Plot the original time series
    #plt.plot(ts.time.jd, ts['flux'], '.', markersize=1, color='black')

    # Fold the time series with the specified period and reference MJD  
    epoch_time = Time(reference_mjd, format='mjd')  # Create a Time object for the reference MJD
    folded_ts = ts.fold(period=period * u.day, epoch_time=epoch_time)
        #plt.plot(folded_ts.time.jd, folded_ts['flux'], '.', markersize=1, color='red')
    print(f'folded timeseries: {folded_ts}')
    print( f'length of folded timeseries: {len(folded_ts)}' )
    
    plt.plot(folded_ts.time.jd, folded_ts['flux'], '.', markersize=1, color='red')
    # calculate chi2
    #print(float(period))
    chi2 = get_chi2( ts, folded_ts, period)
    '''
    # Plot the folded time series
    stats = calculate_segment_stats(combined_folded_ts, 10)
    #print(stats[0])
    folded_ts.time.format = 'jd'
    plt.xlabel('Folded Time (JD)')
    plt.ylabel('flux [a.u.]')
    plt.title("Original Time Series")
    plt.figure(figsize=(10,6))
    plt.plot(combined_folded_ts.time.jd, combined_folded_ts['flux'], '.', markersize=1, color='black', label=file_path)
    ts_binned = aggregate_downsample(combined_folded_ts, time_bin_size=0.001 * u.day)
    ts_binned_2 = aggregate_downsample(combined_folded_ts, time_bin_size=0.001 * u.day, aggregate_func=np.std)
    #print(f'mean folded signal: {ts_binned}')
    #print(f'std of folded signal: {ts_binned_2}')
    plt.plot(ts_binned.time_bin_start.jd, ts_binned['flux'], 'r-', drawstyle='steps-post')
    plt.plot(ts_binned_2.time_bin_start.jd, ts_binned['flux'] + ts_binned_2['flux'], '-', color='orange', drawstyle='steps-post')
    plt.plot(ts_binned_2.time_bin_start.jd, ts_binned['flux'] - ts_binned_2['flux'], '-', color='orange', drawstyle='steps-post')
    #plt.plot(stats[0], stats[1], linestyle='None',marker='.', color='red')
    '''
    plt.xlabel('Folded Time (JD)')
    plt.ylabel('flux [a.u.]')
    
    #plt.legend(loc='upper right')
    plt.title(f'Folded Time Series. Period = {period} days.')
    plt.savefig('folded_timeseries.pdf')
    plt.show()
    
    return chi2

if __name__ == '__main__':
    arguments = docopt(__doc__)
    directory = arguments['--directory']
    period = float(arguments['--period'])
    reference_mjd = float(arguments['--reference_mjd'])
    #combined_ts = combine_timeseries(directory)
    #mean, std = get_stats_timeseries( combined_ts )
    chi2 = analyze_timeseries(directory, period, reference_mjd)
    #chi2 = analyze_split_timeseries(directory, period, reference_mjd)
    print( f'chi2: {chi2}' )









    
# #!/usr/bin/env python3
# """
# Analyze time series data and create a combined plot.

# Usage:
#   analyze_timeseries.py --directory=<directory>
#   analyze_timeseries.py (-h | --help)

# Options:
#   -h --help              Show this help message and exit.
#   --directory=<directory>  Directory where timeseries files are stored.

# """

# import os
# import glob
# import numpy as np
# from astropy.timeseries import TimeSeries
# import matplotlib.pyplot as plt
# from docopt import docopt

# def analyze_timeseries(directory):
#     # Find all files named 'timeseries_*.dat' in the specified directory
#     file_pattern = os.path.join(directory, 'timeseries_*.dat')
#     timeseries_files = glob.glob(file_pattern)

#     if not timeseries_files:
#         print(f"No timeseries files found in the directory: {directory}")
#         return

#     # Sort files by the integer value i in the filename
#     timeseries_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))

#     # Create a combined plot
#     plt.figure(figsize=(10, 6))
#     for file_path in timeseries_files:
#         ts = TimeSeries.read(file_path, format='ascii', time_column='time', time_format='mjd')
#         plt.plot(ts.time.mjd, ts['flux'], '.', markersize=1, label=file_path)

#     plt.xlabel('Julian Date')
#     plt.ylabel('flux')
#     plt.legend(loc='upper right')
#     plt.show()

# if __name__ == '__main__':
#     arguments = docopt(__doc__)
#     directory = arguments['--directory']
#     analyze_timeseries(directory)









# from astropy.timeseries import TimeSeries
# import matplotlib.pyplot as plt
# ts = TimeSeries.read('timeseries.dat', format='ascii', time_column = 'time', time_format = 'mjd')
# plt.plot(ts.time.mjd, ts['flux'], 'k.', markersize = 1)
# plt.xlabel('Julian Date')
# plt.ylabel('flux')
# plt.show()
