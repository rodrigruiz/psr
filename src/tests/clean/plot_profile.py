#!/usr/bin/env python3
"""
Plot a table from an ascii.ecsv file.

Usage:
  script.py --filepath=<filepath> --output=<output> --phasebins=<phasebins> --profile=<profile> [--color=<color>] [--label=<label>] [--title=<title>]
  script.py (-h | --help)

Options:
  -h --help               Show this help message and exit.
"""

from docopt import docopt
#from astropy.timeseries import TimeSeries
from stingray.pulse.search import plot_profile
from add_folded_profiles import read_table
import matplotlib.pyplot as plt

plt.style.use('~/software/psr/src/stingray/clean/latex.mplstyle')

#def plot_pro(file_path, x_columnname, y_columnname, output='timeseries', color='blue', label=None, title=None):
    
if __name__ == '__main__':
    arguments = docopt(__doc__)
    filepath = arguments['--filepath']
    output = arguments['--output']
    phasebins = arguments['--phasebins']
    profile = arguments['--profile']
    
    data = read_table(filepath)
    fig, ax = plt.subplots()
    plot_profile(data[phasebins], data[profile], ax=ax)
    fig.savefig(output + '.png', bbox_inches='tight')
    