""" Chi2s

Usage: Plot.py <INPUT_FILES>... [--wildcard] -o OUTPUT_DIR (--chi2profile | --foldedprofile | --run_hist --runnumber=<runnumber>) [--filepattern=<filepattern>] [--nbin=<int>] [--frequencies=<frequencies>] [--bins_hist=<int>] [--latex] [--signal_strength=<signal_strength>]

Options:
  -h --help                              Help
  <INPUT_FILES>...           Input files
  -o --output_dir OUTPUT_DIR             Output directory
     --nbin=<int>                        Number of bins in the folded profile [default: 32]
     --bins_hist=<int>                   Number of bins in the histogram [default: 200]
     --latex                             Whether to use latex formatting. [default: False]

"""
#  python3 antares_plot_chi2.py -i './folded/Antares_*_chi2_*.hdf5' -o './'
import os, glob
from docopt import docopt
import re
import h5py
import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('~/software/Psr/src/latex.mplstyle')
import plens.TimeSeries as TS
#import scripts 
from scipy.stats import chi2
from epochfolding.stingray_epochfolding import plot_efstat
from epochfolding.stingray_epochfolding import epochfolding_single
from plens.plot_latex_size import set_size
width = 418.25368
figsize = set_size(width, fraction=1, ratio='golden')

def plotHistogram(data, chi2s, output):
    """ Plots a histogram over observed chi2 values and a theoretical chi2 distribution with nbin-1 d.o.f.
    
    Parameters
    ----------
        data : dictionary
            Dictionary containing the command-line arguments.
            To build the correct structure, do:
            data = {}
            for key in arguments:
                data[key.replace("-", "")] = arguments[key] 
        
        chi2s : list or array-lile
            The chi2s values to bin.
        
        output : string
            The filepath for the outputfile.
        
    Returns
    -------
        Saves the plot to `output`.
    
    """
    fig, ax = plt.subplots(figsize=figsize)
    plt.hist(chi2s, bins=int(data['bins_hist']), density=True, alpha=0.6, color='b', label='Observed')
     
    # Plot a chi-squared distribution with appropriate degrees of freedom
    h_axis = np.linspace(0, 200, 100)
    #v_axis = chi2.pdf(h_axis, int(data['nbin']) - 4)
    
    plt.plot(h_axis, chi2.pdf(h_axis, int(data['nbin']) - 1), 'r-', lw=2, label=r'$\chi^2(\mathrm{{d.o.f.}} = {})$'.format(int(data['nbin']) - 1))

    plt.yscale('log')
    plt.xlabel(r'$\chi^2$')
    plt.ylabel('Probability Density')
    plt.xlim([0,200])
    plt.legend()
        
    fig.savefig(data['output_dir'] + 'hist.png',  bbox_inches='tight')
    #plt.show()
    plt.close()
    
    
def getOutputFilepath(data, file=None):
    """ Builds the outputfilepath depending on the command-line arguments chosen.
    
    Parameters
    ----------
        data : dictionary
            Dictionary containing the command-line arguments.
            To build the correct structure, do:
            data = {}
            for key in arguments:
                data[key.replace("-", "")] = arguments[key]
                
    Returns
    -------
        output : string
            The filepath for the outputfile.
            
    """
    
    if data['filepattern']:
        split = re.split(data['filepattern'], file)
        run_number = split[1]
        #gtis['run_number'] = 
    #print(split)
    
        if split[2] != '':
            split_number = split[2]
        #print('Processing Run Nr.: ' + str(run_number) + ', split: ' + str(split_number), end='\n')
            if data['chi2profile']:
                
                out = data['output_dir'] + 'Antares_' + run_number + '_' + split_number + '_chi2profile'
                if data['signal_strength'] is not None:
                    out += '_' + data['signal_strength']
                return out
                
            if data['foldedprofile']:
                return 'Antares_' + run_number + '_' + split_number + '_foldedprofile'
            
        elif data['chi2profile']:
            out =  data['output_dir'] + 'Antares_' + run_number + '_chi2profile'
            if data['signal_strength'] is not None:
                out += '_' + data['signal_strength']
            return out
        
        elif data['foldedprofile']:
                return 'Antares_' + run_number + '_foldedprofile'
            
    elif data['chi2profile']:
        out = data['output_dir'] + 'Antares_chi2profile'
        if data['signal_strength'] is not None:
            out += '_' + data['signal_strength']
        return out
    
    elif data['foldedprofile']:
        return 'Antares_foldedprofile'

def plotRunHist(arguments, data_to_bin, bins, runnumber=0):
    '''the data_to_bin has to be in Hz'''
    fig, ax = plt.subplots(figsize=figsize)
    plt.hist(data_to_bin/1000., bins=bins, density=True, alpha=0.6, color='b', label='ANTARES rates'+ '\n' + 'run number {}'.format(runnumber)+ '\n'+ 'run length {}\,s'.format(round(len(data_to_bin)*0.1, 2)))
    plt.axvline(np.mean(data_to_bin/1000.), color='#d62728', ls='--', label='Mean rate ${}\,\mathrm{{kHz}}$'.format(round(np.mean(data_to_bin/1000.))))
    plt.xlabel(r'ANTARES rates [kHz]')
    plt.ylabel('Probability Density')
    #plt.xlim([0,200])
    plt.legend()
        
    fig.savefig(arguments['output_dir'] + 'run_hist_' + runnumber + '.pdf',  bbox_inches='tight')
    #plt.show()
    plt.close()
    
    print(len(data_to_bin)*np.mean(data_to_bin))
    
def main():
    arguments = docopt(__doc__)

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
   
    input_files.sort()

    for file in input_files:
        output = getOutputFilepath(data, file)
            
        if data['chi2profile']:
            with h5py.File(file) as f: 
                plot_efstat(f['efstats/frequencies'][()], f['efstats/chi2s'][()], nbin=int(data['nbin']), output=output, label=r'$\chi^2$ landscape', true_frequency=None)
                
        if data['foldedprofile']:
            # input is eventlist, not the saved folded profile
            frequencies = list( map(float, data['frequencies'].strip('[]').split(',')) )
            epochfolding_single(file, frequencies, nbin=int(data['nbin']), output=output, plot=True, save=False, format='hdf5', outputdir=data['output_dir'])
        
        if data['run_hist']:
            with h5py.File(file) as h5_file:
                ts, timeslice_duration = TS.readTimeSeries(h5_file)
                plotRunHist(data, ts[ ts['rateOn'] > 0.]['rateOn'].value, int(data['bins_hist']), runnumber=data['runnumber'])
                

if __name__ == "__main__":
    main()