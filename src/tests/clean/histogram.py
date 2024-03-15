"""
Plot Chi2 pdf.
"""

import os, glob
import numpy as np
from docopt import docopt
from astropy.table import Table
from astropy.io import ascii
import matplotlib.pyplot as plt
from add_folded_profiles import read_table

data = read_table('chi2_added_.dat')
fig, ax = plt.subplots()
_ = plt.hist(data['chi2'], bins=20, density=True, color='blue')
#hist, bin_edges = np.histogram(data['chi2'], bins=np.arange(5), density=True)
#print(bin_edges)
#plt.plot(bin_edges, hist)
fig.savefig('chi2_hist.png', bbox_inches='tight')