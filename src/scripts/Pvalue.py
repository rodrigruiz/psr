"""P-value.
Usage: Pvalue.py <INPUT_FILES>... [--wildcard] -o OUTPUT_DIR --filepattern=<filepattern> [--latex] [--nbin=<nbin>]

Options:
    -h help
    -i --input_files INPUT_FILES           Input files
    -o --output_dir OUTPUT_DIR             Output directory
    --nbin=<nbin>                          Number of bins in the folded profile [default: 32]
"""
# python3 Pvalue.py 'Antares_chi2_added*.hdf5' --wildcard -o ./ --filepattern 'Antares_chi2_added(\d*\.*\d+).hdf5'
# the complicated regular expression find integers and floats
import os, glob
from docopt import docopt
import re
import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from scipy.stats.distributions import chi2
from plens.plot_latex_size import set_size
width = 418.25368
figsize = set_size(width, fraction=0.8, ratio='golden')


def main():
    arguments = docopt(__doc__)
    
    #print(stingray.__version__)
    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    print(data)
    
    if data['latex']:
        plt.style.use('~/software/Psr/src/latex.mplstyle')
        
    if data['wildcard']:
        input_files = glob.glob(data['<INPUT_FILES>'][0])
    else:
        input_files = data['<INPUT_FILES>']
        #gti_files = data['gti_files']
    input_files.sort()
    print(input_files)
    
    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])
    
    signal_strength = []
    max_chi2s = []
    max_f = []
    p_values = []
    for file in input_files:
        split = re.split(data['filepattern'], file)
        print(split)
        signal_strength.append(float(split[1]))
        
        with h5py.File(file) as f:
            #print(f['efstats/frequencies'][()], f['efstats/chi2s'][()])
            max_chi2s.append( np.max( f['efstats/chi2s'][()] ) ) 
            max_f.append( f['efstats/frequencies'][()][np.where(f['efstats/chi2s'][()] == np.max( f['efstats/chi2s'][()] ) )][0] )
            p_values.append( chi2.sf(np.max( f['efstats/chi2s'][()] ) , int(data['nbin']) -1) )
    print(signal_strength, max_chi2s, max_f, p_values)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    plt.plot(signal_strength, p_values, 'o')
    plt.axhline(y=2.87e-7, color='r', alpha=0.5, label=r'5$\sigma$' )
    plt.xlabel(r'Signal strength $\mu$ [a.u.]')
    plt.ylabel(r'$p$-value')
    #plt.xlim([0,200])
    plt.legend()
    #ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.set_ylim(np.min(p_values), np.max(p_values))
    #"""
    #ax2 = ax.twinx()
    #alignYaxes([ax,ax2])#,[p_values[0],max_chi2s[-1]])
    #ax2.set_ylim(np.max(max_chi2s), np.min(max_chi2s))
    
    indices = [0, 5,6,7,9]
    max_chi2s_filter = [max_chi2s[x] for x in indices]
    p_values_filter = [p_values[x] for x in indices]
    max_chi2s_strings = [str(round(i,1)) for i in max_chi2s_filter]
    #ax2.set_yticklabels(max_chi2s_strings)
    #ax.get_shared_y_axes().join(ax, ax2)
    ax2 = ax.secondary_yaxis('right')
    ax2.set_yticks(p_values_filter, max_chi2s_strings)
    #ax2.set_yticklabels(max_chi2s_strings)
    #ax2.set_yscale('log')
    ax2.set_ylabel(r'$\chi^2_{\mu,\mathrm{obs}}$')
    #"""
    
    """
    axins = zoomed_inset_axes(ax, 3, loc=1)
    x1, x2, y1, y2 = 0.068, 0.102, -3e-3, 3e-3
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    axins.yaxis.get_major_locator().set_params(nbins=2)
    axins.xaxis.get_major_locator().set_params(nbins=4)
    axins.tick_params(labelleft=False, labelbottom=False)
    axins.plot(signal_strength, p_values, '+')
    axins.axhline(y=2.87e-7, color='r', alpha=0.5 )
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    """
    fig.savefig(data['output_dir'] + 'pvalues.pdf',  bbox_inches='tight')
    #plt.show()
    plt.close()
    
    
if __name__ == "__main__":
    main()
