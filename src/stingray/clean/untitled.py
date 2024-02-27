"""
Plot Chi2 pdf.
"""

import os, glob
import numpy as np
from docopt import docopt
from astropy.table import Table
from astropy.io import ascii
from add_folded_profiles import read_table

data = read_table('chi2_added_.dat')
_ = pl