import os, glob
import numpy as np
from docopt import docopt
from astropy.table import Table
from astropy.io import ascii
from add_folded_profiles import read_table
from analyse_timeseries import plot_chi2_landscape

data = read_table('chi2_added_.dat')
#print(np.max(data['chi2']), data['chi2'].argmax(), data['frequency'][data['chi2'].argmax()] )
plot_chi2_landscape(chi2s=data['chi2'], periods=data['frequency'])