""" Plot Chi2 Distribution.

Usage: antares_plot_chi2.py -i INPUT_FILES -o OUTPUT_DIR [--nbin=<int>] [--bins_hist=<int>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output directory
     --nbin=<int>                        Number of bins in the folded profile [default: 32]
     --bins_hist=<int>                   Number of bins in the histogram [default: 200]

"""
#  python3 antares_plot_chi2.py -i './folded/Antares_*_chi2_*.hdf5' -o './'
import os, glob
from docopt import docopt
import h5py
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('/home/hpc/capn/mppi148h/software/Psr/src/latex.mplstyle')
import plens.TimeSeries as TS
#import scripts 
from scipy.stats import chi2

def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    input_files = glob.glob(data['input_files'])
    input_files.sort()
    
    chi2s = None
    for file in input_files:
        
        with h5py.File(file) as input_file:
            if chi2s is None:
                chi2s = list(input_file['efstats/chi2s'][()])
            else:
                chi2s += list(input_file['efstats/chi2s'][()])
    
    fig, ax = plt.subplots()
    plt.hist(chi2s, bins=int(data['bins_hist']), density=True, alpha=0.6, color='b', label='Observed')
    
     
        # Plot a chi-squared distribution with appropriate degrees of freedom
    h_axis = np.linspace(0, 200, 100)
    v_axis = chi2.pdf(h_axis, int(data['nbin']) - 4)

    #plt.plot(h_axis, v_axis, 'r-', lw=2, label=r'$\chi^2(\mathrm{{d.o.f.}} = {})$'.format(int(data['nbin']) - 4))
    plt.plot(h_axis, chi2.pdf(h_axis, int(data['nbin']) - 1), 'r-', lw=2, label=r'$\chi^2(\mathrm{{d.o.f.}} = {})$'.format(int(data['nbin']) - 1))
    #b= 40
    #plt.plot(h_axis, chi2.pdf(h_axis, int(b)), 'r--', alpha=0.6, label=r'$\chi^2(\mathrm{{d.o.f.}} = {})$'.format(b))
    plt.yscale('log')
    plt.xlabel('$\chi^2$')
    plt.ylabel('Probability Density')
    plt.xlim([0,200])
    plt.legend()
        
    fig.savefig(data['output_dir'] + 'hist.png',  bbox_inches='tight')
    #plt.show()
    plt.close()


if __name__ == "__main__":
    main()