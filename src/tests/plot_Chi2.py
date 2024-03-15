import os, glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('~/software/psr/src/stingray/clean/latex.mplstyle')
import plens.TimeSeries as TS
#import scripts 
import sys
sys.path.insert(0, '../clean')
#from stingray_epochfolding import epochfolding_scan, plot_efstat
from stingray.events import EventList
from scipy.stats import chi2
chi2s = []
nbin = 31

file = 'chi2s.hdf5'
with h5py.File(file) as f:
    chi2s = list(f['efstats/chi2s'][()])
    
print(len(chi2s))
#print(chi2s)
fig, ax = plt.subplots()
plt.hist(chi2s, bins=200, density=True, alpha=0.6, color='b', label='Observed')

     
        # Plot a chi-squared distribution with appropriate degrees of freedom
h_axis = np.linspace(0, 200, 100)
v_axis = chi2.pdf(h_axis, nbin)
plt.plot(h_axis, v_axis, 'r-', lw=2, label=f'Chi-squared df={nbin}')
plt.yscale('log')
plt.xlabel('Chi-squared')
plt.ylabel('Probability Density')
plt.xlim([0,200])
plt.legend()
        
fig.savefig('hist.png',  bbox_inches='tight')
#plt.show()
plt.close()
#"""