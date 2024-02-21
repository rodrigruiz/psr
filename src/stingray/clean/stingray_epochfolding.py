#!/usr/bin/env python3
"""
Epochfolding.

Usage:
  stingray_epochfolding.py  -i INPUT_FILE --obs_length=<obs_length> --testperiods=<testperiods> [--mode=<mode>] [--nbin=<nbin>] [--oversampling=<oversampling>] [--true_period=<true_period>] [--plot=<plot>] [--save=<save>] [--root_dir=<root_dir>] [--fit=<fit>] [--pulseshape=<pulseshape> ] 
  stingray_epochfolding.py (-h | --help)

Options:
  -h --help                   Show this help message and exit.
  -i INPUT_FILE               Filepath to the eventlist (in the ascii.escv format).
  --mode=<mode>               Which functions to run: 'single' (then `testperiods` should be set),
                                                      'scan' (then either `testperiods` or `true_period` should be set),
                                                      'all' (both `testperiods` and `true_period` should be set), [default: all]
  --nbin=<nbin>               Number of bins in the folded timeseries. [default: 32]
  --oversampling=<oversampling>  Oversampling of the folded timeseries. [default: 10]
  --obs_length=<obs_length>   Observation length.
  --testperiods=<testperiods> Testperiods to perform the epoch folding on. Can also be only one value or `true_period`
  --true_period=<true_period> True period of the signal (if known).
  --plot=<plot>               Whether to plot the results [default: True]
  --save=<save>               Whether to save **currently** the output of epochfolding_single()
  --root_dir=<root_dir>       root directory where to save the folded pulse profiles
  --fit=<fit>                 Whether to fit the epoch folding statistics [default: True]
  --pulseshape=<pulseshape>   Pulseshape to fit the statistics to. 
                              'arbitrary' will fit a gaussian, 
                              'sine' will fit a sinc function [default: arbitrary]
  
"""
import os
from docopt import docopt
from astropy.table import Table
from astropy.io import ascii
from collections.abc import Iterable
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('../../latex.mplstyle')
import seaborn as sb
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (10, 6)

#from generate_split_timeseries import generate_pulse_train_gauss, gauss
#from stingray import Lightcurve
from stingray.events import EventList
from stingray.pulse.search import epoch_folding_search
from stingray.pulse.pulsar import fold_events
from stingray.stats import fold_detection_level
from stingray.pulse.search import plot_profile, search_best_peaks
from stingray.pulse.modeling import fit_gaussian, fit_sinc

def epochfolding_single(filepath, testperiods, nbin=32, output='profile', plot=True, save=True):
    """Performs epochfolding with a certain testperiod or on a range of testperiods.
     
    Parameters
    ----------
        filepath : str 
            Filename or entire path of the eventlist.
        
        testperiods : float or list of floats
            Testperiod(s) to test with the epochfolding algorithm.
        
        nbin : int
            Number of bins in the folded eventlist.
            Defaults to 32.
          
    Optional Parameters
    -------------------
        plot : bool
            Whether to plot the folded pulse profiles.
            
        output : str
            filename of the plot.
            
    Returns
    -------
        phase_bins : array of floats
            The phases corresponding to the pulse profile

        profile : array of floats
            The pulse profile (folded counts)

        profile_err : array of floats
            The uncertainties on the pulse profile
        
    """
    
    events = EventList().read(filepath, 'ascii')
    
    fig, ax = plt.subplots()
    
    # perform epochfolding and plot resulting profile
    if isinstance(testperiods, Iterable):
        for p in testperiods:
            phase_bins, profile, profile_err = _fold_events(events.time, p, nbin=nbin, save=save)
            if plot:
                ax = plot_profile(phase_bins, profile, ax=ax)
                fig.savefig(output + '.png', bbox_inches='tight')
            '''
            phase_bins, profile, profile_err = fold_events(events.time, 1/p, nbin=nbin)
            if plot:
                ax = plot_profile(phase_bins, profile, ax=ax)
                fig.savefig(output + '.png', bbox_inches='tight')
            if save:
                save_to_ascii(phase_bins, profile, profile_err, p)
                
                # defining astopy table
                data = Table()
                data['phase_bins'] = phase_bins
                data['folded_profile'] = profile
                data['folded_profile_err'] = profile_err
                ascii.write(data, 'folded_profile_{}.dat'.format(p), overwrite=True)
                
            '''
    elif isinstance(testperiods, float):
        phase_bins, profile, profile_err = _fold_events(events.time, testperiods, nbin=nbin, save=save)
        if plot:
            ax = plot_profile(phase_bins, profile, ax=ax)
            fig.savefig(output + '.png', bbox_inches='tight')
        '''
        phase_bins, profile, profile_err = fold_events(events.time, 1/testperiods, nbin=nbin)
        if plot:
            ax = plot_profile(phase_bins, profile, ax=ax)
            fig.savefig(output + '.png', bbox_inches='tight')
        if save:
            save_to_ascii(phase_bins, profile, profile_err, testperiods)
        '''
    else:
        raise ValueError(
            "testperiods should be either a float of a list of floats!"
        )
    
        
    return phase_bins, profile, profile_err

def _fold_events(times, testperiod, nbin=32, save=False, root_dir=None):
    """Helping function.
       Calculate the folded pulse profile.
    """
    
    phase_bins, profile, profile_err = fold_events(times, 1/testperiod, nbin=nbin)
    if save:
        save_to_ascii(phase_bins, profile, profile_err, testperiod, root_dir=root_dir)
    return phase_bins, profile, profile_err

def save_to_ascii(phase_bins, profile, profile_err, testperiod, root_dir=None):
    """Helping function.
       Saves the output of fold_events() to an ascii file.
    """
    if root_dir is None:
        output_dir = './folded_events/'
    else:
        output_dir = './folded_events/' + root_dir
    
    if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        
    data = Table()
    data['phase_bins'] = phase_bins
    data['folded_profile'] = profile
    data['folded_profile_err'] = profile_err
    ascii.write(data, output_dir + 'folded_profile_{}.dat'.format(testperiod), overwrite=True)
    
def epochfolding_scan(filepath, testperiods, nbin=32, oversampling=10, obs_length=None, plot=True, true_period=None):
    """Performs epochfolding for a range of testperiods.
    
    Parameters
    ----------
        filepath : str
            Filename or entire path of the eventlist.
            
        testperiods : float or list of floats
            Testperiods to test with the epochfolding algorithm.
            If a single float, a narrow range around the given period will be performed.
            For epochfolding with a single period, please use epochfolding_single().
            
    Optional Parameters
    -------------------
        nbin : int
            Number of bins in the folded eventlist.
            Defaults to 32.
            
        oversampling : int
            Oversampling of the folded eventlist.
            Defaults to 10.
            To be specified when `testperiods` is a single float.
            Together with `obs_length` determines the precision of the testperiods.
        
        obs_length : float
            Observation length
            To be specified when `testperiods` is a single float.
            Together with `oversampling` determines the precision of the testperiods.
        
        plot : bool
            Whether to plot the statistics of the epoch folding search.
            
        true_period : float
            True period of the signal (if known).
            
    Returns
    -------
        frequencies : array-like
            epoch folding trial periods.
            
        effreq : array-like
            frequency grid of the folded timeseries.
            
        efstat : array-like
            epoch folding statistics (chi2-values) for the frequency bins in `effreq`.
    
    """
    events = EventList().read(filepath, 'ascii')
        
    if isinstance(testperiods, float): 
        if obs_length is None:
            raise ValueError(
                "Please specify the observation length if you enter only one testperiod."
            )
            
        df_min = 1/obs_length
        df = df_min / oversampling
        frequencies = np.arange(1/testperiods - 200 * df, 1/testperiods + 200 * df, df)
    elif isinstance(testperiods, Iterable):
        frequencies = 1/testperiods
    else:
        raise ValueError(
            "testperiods should be a float or a list of floats!"
        )

    effreq, efstat = epoch_folding_search(events.time, frequencies, nbin=nbin)
    
    if plot:
        plot_efstat(effreq, efstat, nbin, true_period=true_period)
    
    return frequencies, effreq, efstat


def plot_efstat(effreq, efstat, nbin, output='chi2', label='EF statistics', true_period=None):
    """Plots the chi2 distribution of the epochfolding scan in frequency space.
    
    Parameters
    ----------
        effreq : array-like
            frequency grid of the folded timeseries.
            
        efstat : array-like
            epoch folding statistics (chi2-values) for the frequency bins in `effreq`.
        
        nbin : int
            Number of bins in the folded eventlist.
            
    Optional Parameters
    -------------------
        
        output : str
            Output filename or entire path to the outputfile (without the fileending!).
            Defaults to 'chi2'.

        label : str
            Label to print in legend.
            Defaults to 'EF statistics'.
            
        true_period : float
            True period of the signal (if known).
            Default is None.
        
    """
    
    plt.figure()
    plt.plot(effreq, efstat, label=label, color='blue')
    plt.axhline(nbin - 1, ls='--', lw=2, color='k', label='n - 1')
    if isinstance(true_period, float):
        plt.axvline(1/true_period, lw=3, alpha=0.5, color='r', label='Correct frequency $f_\mathrm{true} = {' + str(1/true_period) + '}$')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('EF Statistics')
    plt.legend()
    plt.savefig(output + '.png', bbox_inches='tight')
    

def fit_stats( frequencies, nbin, effreq, efstat, obs_length, pulseshape='arbitrary', plot=False, true_period=None ):
    """Fits a gaussian or a sinc-function to the epoch folding statistics around ***the first peak***.
    
    Parameters
    ----------
        frequencies : array-like
            epoch folding trial periods.
            
        nbin : int
            Number of bins in the folded eventlist.
            
        effreq : array-like
            frequency grid of the folded timeseries.
            
        efstat : array-like
            epoch folding statistics (chi2-values) for the frequency bins in `effreq`.
            
        obs_length : float
            Observation length.
            
        pulseshape : str
            if 'arbitrary' a gauss function is fit to the epoch folding statistics (default).
            if 'sine' a sinc function is fit.
            
    Optional Parameters
    -------------------
        plot : bool
            Whether to plot the statistics of the epoch folding search.
            
        true_period : float
            True period of the signal (if known).
            
    Returns
    -------
        cand_freq_ef : array-like
            array containing the x position of the peaks above threshold. If no
            peaks are above threshold, an empty list is returned. The array is
            sorted by inverse value of statistics.
            
        cand_stat_ef :
            for each peak in `cand_freq_ef`, give the corresponding chi2 value. Empty if no peaks
            above threshold.
            
        fit : function
            The best-fit function, accepting x as input
            and returning the best-fit model as output.
        
    """
    
    # search best candidate frequencies
    df_min = 1/obs_length
    ntrial = (frequencies[-1] - frequencies[0]) / df_min
    ef_detlev = fold_detection_level(nbin, epsilon=0.001, ntrial=len(effreq))
    cand_freqs_ef, cand_stat_ef = search_best_peaks(effreq, efstat, ef_detlev)
    
    if pulseshape == 'arbitrary':
        fit = fit_gaussian(effreq, efstat-(nbin-1),amplitude=max(efstat-(nbin-1)), 
                mean=cand_freqs_ef[0], stddev=1/(np.pi*obs_length))
    elif pulseshape == 'sine':
        fit = fit_sinc(effreq, efstat-(nbin-1),amp=max(efstat-(nbin-1)), mean=cand_freqs_ef[0], 
            obs_length=obs_length)
    else:
        raise ValueError(
            "pulseshape should be either `arbitrary` or `sine`!"
        )
        
    if plot:
        plot_stats_fit(fit, effreq, efstat, nbin, output='chi2_fit', label='EF statistics', true_period=true_period)
    
    return cand_freqs_ef, cand_stat_ef, fit

def plot_stats_fit(fit, effreq, efstat, nbin, output='chi2_fit', label='EF statistics', true_period=None):
    """Plots the chi2 distribution of the epochfolding scan in frequency space.
    
    Parameters
    ----------
        fit : function
            The best-fit function, accepting x as input
            and returning the best-fit model as output.
        
        effreq : array-like
            frequency grid of the folded timeseries.
            
        efstat : array-like
            epoch folding statistics (chi2-values) for the frequency bins in `effreq`.
        
        nbin : int
            Number of bins in the folded eventlist.
            
    Optional Parameters
    -------------------
        
        output : str
            Output filename or entire path to the outputfile (without the fileending!).
            Defaults to 'chi2_fit'.

        label : str
            Label to print in legend.
            Defaults to 'EF statistics'.
            
        true_period : float
            True period of the signal (if known).
            Default is None.
        
    """
        
    plt.figure(figsize=(15, 5))
    plt.plot(effreq, efstat-(nbin-1), label=label, color='blue')
    plt.plot(effreq, fit(effreq), label='Best fit', color='red', ls='--')
    if isinstance(true_period, float):
        plt.axvline(1/true_period, alpha=0.5, color='r', label='Correct frequency $f_\mathrm{true} = {' + str(1/true_period) + '}$')
    plt.axvline(fit.mean[0], alpha=0.5, label=r'Fit frequency $f_\mathrm{fit} = $' + '{}'.format(fit.mean[0]))
    #plt.xlim( [fg.mean[0] * 0.95, fg.mean[0] * 1.05] )
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('EF Statistics')
    plt.legend()
    plt.savefig(output + '.png', bbox_inches='tight')
    


'''
#residuals
#plt.figure(figsize=(15, 5))
plt.plot(freq, efstat-(nbin-1)-fg(freq))
plt.xlabel('Frequency (Hz)')
_ = plt.ylabel('Residuals')

'''

if __name__ == '__main__':
    arguments = docopt(__doc__)
    filepath = arguments['-i']
    mode = arguments['--mode']
    nbin = int(arguments['--nbin'])
    #oversampling = int(arguments['--oversampling'])
    obs_length = float(arguments['--obs_length'])
    if arguments['--testperiods'][0] != '[':
        testperiods = float(arguments['--testperiods'])
    else: 
        testperiods = map(float, arguments['--testperiods'].strip('[]').split(','))
    #true_period = float(arguments['--true_period'])
    plot = bool(arguments['--plot'])
    save = bool(arguments['--save'])
    #fit = bool(arguments['--fit'])
    #pulseshape = arguments['--pulseshape']
    
    if mode == 'single':
        if nbin is None:
             epochfolding_single(filepath, testperiods, output='profile', plot=plot, save=save)
        else:
            epochfolding_single(filepath, testperiods, nbin=nbin, output='profile', plot=plot, save=save)
        
    elif mode == 'scan':
        frequencies, effreq, efstat = epochfolding_scan(filepath, true_period, nbin=nbin, oversampling=oversampling, 
                          obs_length=obs_length, plot=plot, true_period=true_period)
        if fit:
            fit_stats( frequencies, nbin, effreq, efstat, obs_length, pulseshape=pulseshape, plot=plot, true_period=true_period )
            
    elif mode == 'all':
        epochfolding_single(filepath, testperiods, nbin=nbin, output='profile', plot=plot, save=save)
        frequencies, effreq, efstat = epochfolding_scan(filepath, true_period, nbin=nbin, oversampling=oversampling,
                                                        obs_length=obs_length, plot=plot, true_period=true_period)
        cand_freqs_ef, cand_stat_ef, fit = fit_stats( frequencies, nbin, effreq, efstat, obs_length,
                                                     pulseshape=pulseshape, plot=plot, true_period=true_period )
        