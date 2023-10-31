#!/usr/bin/env python3
"""
Generate and split time series into multiple files.

Usage:
  script.py --time_start=<time_start> --time_stop=<time_stop> --number_of_points=<number_of_points> --signal_strength=<signal_strength> --period=<period> --num_files=<num_files>
  script.py (-h | --help)

Options:
  -h --help               Show this help message and exit.
  --time_start=<time_start>     Start time
  --time_stop=<time_stop>       Stop time
  --number_of_points=<number_of_points> Number of points
  --signal_strength=<signal_strength>   Signal strength
  --period=<period>           Period of the sinusoid signal
  --num_files=<num_files>     Number of output files
"""

import numpy as np
from astropy.timeseries import TimeSeries
from astropy.time import Time
from docopt import docopt

def generate_time_series(time_start, time_stop, number_of_points, signal_strength, period, num_files):
    times = Time(np.linspace(float(time_start), float(time_stop), int(number_of_points)), format='mjd')
    periodic_signal = float(signal_strength) * np.sin(2 * np.pi * times.jd / float(period))
    noise = np.random.normal(0, 0.1, int(number_of_points))
    flux = periodic_signal + noise

    # Split the time series into multiple files
    file_indices = np.array_split(np.arange(int(number_of_points)), int(num_files))
    for i, indices in enumerate(file_indices):
        ts = TimeSeries(time=times[indices])
        ts['flux'] = flux[indices]
        file_name = f'timeseries_{i}.dat'
        ts.write(file_name, format='ascii')

if __name__ == '__main__':
    arguments = docopt(__doc__)
    generate_time_series(
        arguments['--time_start'],
        arguments['--time_stop'],
        arguments['--number_of_points'],
        arguments['--signal_strength'],
        arguments['--period'],
        arguments['--num_files']
    )
