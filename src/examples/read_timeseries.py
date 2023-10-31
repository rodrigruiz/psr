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


def analyze_timeseries(directory, period, reference_mjd):
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

        # Fold the time series with the specified period and reference MJD
        
        epoch_time = Time(reference_mjd, format='mjd')  # Create a Time object for the reference MJD
        folded_ts = ts.fold(period=period * u.day, epoch_time=epoch_time)
    
        # Sum the folded time series to the combined_folded_ts
        if combined_folded_ts is None:
            combined_folded_ts = folded_ts
        else:
            combined_folded_ts['flux'] += folded_ts['flux']

    # Plot the folded time series
    folded_ts.time.format = 'jd'
    plt.plot(combined_folded_ts.time.jd, combined_folded_ts['flux'], '.', markersize=1, label=file_path)
    plt.xlabel('Folded Time (MJD)')
    plt.ylabel('flux')
    plt.legend(loc='upper right')
    plt.title(f'Combined Folded Time Series')
    plt.show()

if __name__ == '__main__':
    arguments = docopt(__doc__)
    directory = arguments['--directory']
    period = float(arguments['--period'])
    reference_mjd = float(arguments['--reference_mjd'])
    analyze_timeseries(directory, period, reference_mjd)









    
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
