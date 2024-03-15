#!/usr/bin/env python3
"""
Plot a timeseries from an ascii.ecsv file.

Usage:
  script.py --file_path=<file_path> [--output=<outputfilepath>] [--data_column_name=<data_column_name>] [--color=<color>] [--label=<label>] [--title=<title>]
  script.py (-h | --help)

Options:
  -h --help               Show this help message and exit.
  --file_path=<file_path> Filepath to the timeseries file.
"""

from docopt import docopt
from astropy.timeseries import TimeSeries
import matplotlib.pyplot as plt
plt.style.use('~/software/psr/src/stingray/clean/latex.mplstyle')

def plot_timeseries(file_path, output='timeseries', time_column_name='time', data_column_name='counts', color='blue', label=None, title=None):
    """Plots a astropy.timeseries.TimeSeries object.
    
    Parameters
    ----------
        file_path : str
            Path to the table file containing the astropy.timeseries.TimeSeries.
    
    Optional Parameters
    -------------------
        
        output : str
            Output filename or entire path to the outputfile (without the fileending!).
            Defaults to 'timeseries'.
        
        data_column_name : str
            Column name of the data to plot.
            Defaults to 'counts'.
        
        color : str
            Color of the data in the plot.
            Defaults to 'blue'.
        
        label : str
            Label to print in legend.
            Defaults to 'Timeseries'.
        
        title : str
            Title of the figure.
            Defaults to 'Timeseries'.
        
    """
    
    ts = TimeSeries.read(file_path, format='ascii', time_column='time', time_format='mjd')
    
    if label is None:
        label = 'Timeseries'
    if title is None:
        title = 'Timeseries'
        
    plt.figure(figsize=(10, 5))
    plt.plot(ts.time.value, ts[data_column_name], label=label, color=color )
    plt.xlabel('Time [MJD]')
    plt.ylabel('Counts [a.u.]')
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.savefig( output + '.png', bbox_inches='tight' )
    plt.show()
    
if __name__ == '__main__':
    arguments = docopt(__doc__)
    plot_timeseries(
        arguments['--file_path'],
        arguments['--output']
    )