#!/usr/bin/env python
"""
Reads ANTARES efficiency files from a directory, and creates an astropy QTable with the coumns:
 - unix timestamp
 - year
 - month
 - day
 - mean efficiency
 - median efficiency
 - standard deviation efficiency
 - standard error efficiency
 - number of PMTs with 0 efficiency
 - number of PMTs with != 0 efficiency

Usage:
  calculate_efficiencies.py -I <input_directory> -o <output_file>

Options:
  -h --help     Show this help message and exit.
  -I <input_directory> Path to the efficiency files
  -o <output_file> Path and name of the output file

Examples:
  calculate_efficiencies.py /path/to/directory

"""
import os
import re
import docopt
import numpy as np
from datetime import datetime
from astropy.table import QTable
from astropy.table import Table
from astropy.io import fits

def parse_file_name(filename):
    """
    Parse the given filename and extract information about the month, day, and year.

    Parameters:
    - filename (str): Name of the file to parse.

    Returns:
    - dict: Dictionary containing 'month', 'day', and 'year'.
    """
    months = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12}
    pattern = r"([A-Za-z]+)(\d+)_(\d+)\.txt"
    d = {'day': 0, 'month': 0, 'year': 0}
    # Use re.match to extract the month abbreviation, x, and year
    match = re.match(pattern, filename)
    if match:
        d['month'] = months[match.group(1)]
        d['day'] = 1 + 6*(int(match.group(2)) - 1)
        d['year'] = int(match.group(3))
        #print(f"file: {filename}. Month: {d['month']}, day: {d['day']}, Year: {d['year']}")
    else:
        print(f"Error: Filename {filename} does not match the expected format.")
    return d

def get_unix_time(year, month, day, hour = 0, minute = 0, second = 0):
    """
    Calculate the Unix timestamp for a given date and time.

    Parameters:
    - year (int): Year.
    - month (int): Month.
    - day (int): Day.
    - hour (int): Hour (default is 0).
    - minute (int): Minute (default is 0).
    - second (int): Second (default is 0).

    Returns:
    - int: Calculated Unix timestamp.
    """
    input_date = datetime(year, month, day, hour, minute, second)
    return int((input_date - datetime(1970, 1, 1)).total_seconds())

def calculate_efficiency_statistics(filename):
    """
    Calculate statistics (mean, median, std, sterr, zeros, nonzeros) from the last column of the given file.

    Parameters:
    - filename (str): Path to the file.

    Returns:
    - dict: Dictionary containing calculated statistics.
    """
    with open(filename, 'r') as file:
        values = [float(line.split()[-1][:-1]) for line in file]
    # Convert the list of values to a numpy array
    values = [0 if np.isnan(x) else x for x in values]
    values_array = np.array(values)
    mean = np.mean(values_array)
    median = np.median(values_array)
    std = np.std(values_array)
    sterr = std/np.sqrt(len(values_array))
    zeros = np.sum(values_array == 0)
    nonzeros = np.sum(values_array != 0)
    return {'mean': np.mean(values_array), 'median': np.median(values_array), 'std':np.std(values_array), 'sterr': sterr, 'zeros': zeros, 'nonzeros': nonzeros}

def main():
    args = docopt.docopt(__doc__)
    directory = args['<input_directory>']
    output_filename = args['output_file']

    if not os.path.isdir(directory):
        print(f"Error: '{directory}' is not a valid directory.")
        return

    files = os.listdir(directory)

    times = mean = median = std = sterr = zeros = nonzeros = names = year = month = day = np.array([])
    for filename in files:
        if filename.endswith('.txt'):
            try:
                d = parse_file_name(filename)
                if d['year'] != 0:
                    unix_time = get_unix_time(d['year'], d['month'], d['day'])
                    stats = calculate_efficiency_statistics(os.path.join(directory,filename))
                    year = np.append(d['year'], year)
                    month = np.append(d['month'], month)
                    day = np.append(d['day'], day)
                    times = np.append(int(unix_time), times)
                    mean = np.append(stats['mean'], mean)
                    median = np.append(stats['median'], median)
                    std = np.append(stats['std'], std)
                    sterr = np.append(stats['sterr'], sterr)
                    zeros = np.append(stats['zeros'], zeros)
                    nonzeros = np.append(stats['nonzeros'], nonzeros)
            except (ValueError, IndexError):
                print(f"Error processing file: {filename}")

    t = Table([times, year, month, day, mean, median, std, sterr, zeros, nonzeros],
               names=['unix_time', 'year', 'month', 'day', 'mean_efficiency', 'median_efficiency', 'std_efficiency', 'sterr_efficiency', 'zero_efficiency', 'non_zero_efficiency'])

    t.sort('unix_time')
    fits.writeto(output_filename, np.array(t), overwrite=True)

if __name__ == '__main__':
    main()
