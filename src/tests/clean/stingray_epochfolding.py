#!/usr/bin/env python3
"""
Epochfolding.

Usage:
  stingray_epochfolding.py  -i INPUT_FILE [--obs_length=<obs_length>] [--frequencies=<frequencies>] [--mode=<mode>] [--nbin=<nbin>] [--oversampling=<oversampling>] [--true_frequency=<true_frequency>] [--plot] [--save] [--outputdir=<outputdir>] [--root_dir=<root_dir>] [--fit --pulseshape=<pulseshape> ] [--number_testf=<number_testf>]
  stingray_epochfolding.py (-h | --help)

Options:
  -h --help                   Show this help message and exit.
  -i INPUT_FILE               Filepath to the eventlist (in the ascii.escv format).
  --mode=<mode>               Which functions to run: 'single' (then `frequencies` should be set),
                                                      'scan' (then either `frequencies` or `true_frequency` should be set),
                                                      'all' (both `frequencies` and `true_frequency` should be set), [default: all]
  --nbin=<nbin>               Number of bins in the folded timeseries. [default: 32]
  --oversampling=<oversampling>  Oversampling of the folded timeseries. [default: 10]
  --obs_length=<obs_length>   Observation length.
  --frequencies=<frequencies> frequencies to perform the epoch folding on. Can also be only one value or `true_frequency`
  --true_frequency=<true_frequency> True frequency of the signal (if known). [default: 10.]
  --plot               Whether to plot the results [default: True]
  --save               Whether to save **currently** the output of epochfolding_single()
  --root_dir=<root_dir>       root directory where to save the folded pulse profiles
  --fit                Whether to fit the epoch folding statistics [default: False]
  --pulseshape=<pulseshape>   Pulseshape to fit the statistics to. 
                              'arbitrary' will fit a gaussian, 
                              'sine' will fit a sinc function [default: arbitrary]
  
"""
import os, glob
import h5py
from docopt import docopt
from astropy.table import Table
from astropy.io import ascii
from collections.abc import Iterable
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('~/software/Psr/src/latex.mplstyle')
import seaborn as sb
import matplotlib as mpl
from plens.plot_latex_size import set_size
width = 418.25368
figsize = set_size(width, fraction=1, ratio='golden')
#mpl.rcParams['figure.figsize'] = set_size(width)

#from generate_split_timeseries import generate_pulse_train_gauss, gauss
#from stingray import Lightcurve
from stingray.events import EventList
from stingray.pulse.search import epoch_folding_search, _folding_search
from stingray.pulse.pulsar import fold_events
from stingray.stats import fold_detection_level
from stingray.pulse.search import search_best_peaks
#from stingray.pulse.search import plot_profile
from stingray.pulse.modeling import fit_gaussian, fit_sinc

def plot_profile(phase, profile, width, err=None, ax=None, fraction=0.5):
    """Plot a pulse profile showing some stats.

    If err is None, the profile is assumed in counts and the Poisson confidence
    level is plotted. Otherwise, err is shown as error bars

    Parameters
    ----------
    phase : array-like
        The bins on the x-axis

    profile : array-like
        The pulsed profile

    Other Parameters
    ----------------
    ax : `matplotlib.pyplot.axis` instance
        Axis to plot to. If None, create a new one.

    Returns
    -------
    ax : `matplotlib.pyplot.axis` instance
        Axis where the profile was plotted.
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=set_size(width, fraction=fraction))

    #mean = np.mean(profile)
    if np.all(phase < 1.5):
        phase = np.concatenate((phase, phase + 1))
        profile = np.concatenate((profile, profile))
    ax.plot(phase, profile, drawstyle="steps-mid")
    #if err is None:
    #    err_low, err_high = poisson_conf_interval(mean, interval="frequentist-confidence", sigma=1)
    #    ax.axhspan(err_low, err_high, alpha=0.5)
    #else:
    #    err = np.concatenate((err, err))
    #    ax.errorbar(phase, profile, yerr=err, fmt="none")
    plt.title('Folded Pulse Profile')
    ax.set_ylabel("Counts [a.u.]")
    ax.set_xlabel(r"Phase $\varphi_\mathrm{ef}$")
    return ax

def load_eventlist(directory):
    file_pattern = os.path.join(directory, 'eventlist_*.dat')
    eventlist_files = glob.glob(file_pattern)
    
    if not eventlist_files:
        print(f"No timeseries files found in the directory: {directory}")
        return
        #FileNotFoundError(f"No timeseries files found in the directory: {directory}")
    
    key = lambda x: int(x.split('_')[-1].split('.')[0])
    eventlist_files.sort(key=key)
    
    return eventlist_files
    '''
    if len(eventlist_files) == 1:
        return EventList().read(filepath, 'ascii')
    else:
        eventlists = []
        for file in eventlist_files:
            Eventlist().read(file, 'ascii')
    '''
        
def epochfolding_single(filepath, frequencies, nbin=32, output='profile', plot=True, save=True, format='ascii', outputdir='./folded_events/'):
    """Performs epochfolding with a certain testfrequency or on a range of testfrequencies.
     
    Parameters
    ----------
        filepath : str 
            Filename or entire path of the eventlist.
        
        frequencies : float or list of floats
            Testfrequency(s) to test with the epochfolding algorithm.
        
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
    if os.path.isfile(filepath):
        events = EventList().read(filepath, format)
        #phase_bins, profile, profile_err = _ef_single(events, testperiods, nbin=nbin, output=output, plot=plot, save=save)
        _ef_single(events, frequencies, nbin=nbin, output=output, plot=plot, save=save, outputdir=outputdir, root_dir='eventlist/')
        
    if os.path.isdir(filepath):
        eventlists = load_eventlist(filepath)
        
        for i, eventlist in enumerate(eventlists):
            events = EventList().read(eventlist, format)
            
            #phase_bins, profile, profile_err = _ef_single(events, testperiods, nbin=nbin, output=output, plot=plot, save=save, root_dir=f'eventlist_{i}/')
            _ef_single(events, frequencies, nbin=nbin, output=output + f'_{i}', plot=plot, save=save, outputdir=outputdir, root_dir=f'eventlist_{i}/')
    
    """   
    fig, ax = plt.subplots()
    
    # perform epochfolding and plot resulting profile
    if isinstance(testperiods, Iterable):
        for p in testperiods:
            phase_bins, profile, profile_err = _fold_events(events.time, p, nbin=nbin, save=save)
            if plot:
                ax = plot_profile(phase_bins, profile, ax=ax)
                fig.savefig(output + '.png', bbox_inches='tight')

    elif isinstance(frequencies, float):
        phase_bins, profile, profile_err = _fold_events(events.time, frequencies, nbin=nbin, save=save)
        if plot:
            ax = plot_profile(phase_bins, profile, ax=ax)
            fig.savefig(output + '.png', bbox_inches='tight')
 
    else:
        raise ValueError(
            "frequencies should be either a float of a list of floats!"
        )
    """
    #return phase_bins, profile, profile_err

def _ef_single(events, frequencies, nbin=32, output='profile', plot=True, save=True, outputdir='./folded_events/', root_dir=None):
    """Helping function.
       Calls _fold_events and plots the folded events (if plot).
    """
    #change fraction here
    fig, ax = plt.subplots(figsize=figsize)
    
    # perform epochfolding and plot resulting profile
    if isinstance(frequencies, Iterable):
        label = []
        for p in frequencies:
            phase_bins, profile, profile_err = _fold_events(events.time, p, nbin=nbin, save=save, outputdir=outputdir, root_dir=root_dir)
            if plot:
                ax = plot_profile(phase_bins, profile, width, ax=ax, fraction=1)
                label.append(r'$f_\mathrm{test} = ' +  str(p) + '\, \mathrm{Hz}$ ')
                ax.legend(label)
                
                fig.savefig(output + '.pdf', bbox_inches='tight')

    elif isinstance(frequencies, float):
        phase_bins, profile, profile_err = _fold_events(events.time, frequencies, nbin=nbin, save=save, outputdir=outputdir, root_dir=root_dir)
        if plot:
            ax = plot_profile(phase_bins, profile, width, ax=ax, fraction=0.5)
            fig.savefig(output + '.pdf', bbox_inches='tight')
    
    else:
        raise ValueError(
            "frequencies should be either a float of a list of floats!"
        )
    
    plt.close()
    print(phase_bins)
    #return phase_bins, profile, profile_err

def _fold_events(times, frequency, nbin=32, save=False, outputdir='./folded_events/', root_dir=None):
    """Helping function.
       Calculate the folded pulse profile.
    """
    
    phase_bins, profile, profile_err = fold_events(times, frequency, nbin=nbin)
    if save:
        save_to_ascii(phase_bins, profile, profile_err, frequency, outputdir=outputdir, root_dir=root_dir)
        
    return phase_bins, profile, profile_err

def save_to_ascii(phase_bins, profile, profile_err, frequency, outputdir='./folded_events/', root_dir=None, output='folded_profile_'):
    """Helping function.
       Saves the output of fold_events() to an ascii file.
    """
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)
        
    if root_dir is None:
        output_dir = outputdir
    else:
        output_dir = outputdir + root_dir

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        
    data = Table()
    data['phase_bins'] = phase_bins
    data['folded_profile'] = profile
    data['folded_profile_err'] = profile_err
    
    #precision = len(str(frequency).split(".")[1])
    #ascii.write(data, output_dir + output + '{0:.{1}f}.dat'.format(frequency, precision), format='ecsv', overwrite=True)
    #print(output_dir + output + '{}.dat'.format(frequency))
    ascii.write(data, output_dir + output + '{0:.4fsave}.dat'.format(frequency), format='ecsv', overwrite=True)

def savehdf5(frequencies, chi2s, filename):
    #filename.create_dataset('efstats/frequencies', frequencies)
    #filename.create_dataset('efstats/chi2s', chi2s)
    filename['efstats/frequencies'] = frequencies
    filename['efstats/chi2s'] = chi2s
    
    return
    
def epochfolding_scan(filepath, frequencies, nbin=32, oversampling=10, number_testf=100, obs_length=None, plot=True, save=False, true_frequency=None, format='ascii', output='chi2s.hdf5'):
    """Performs epochfolding for a range of frequencies.
    
    Parameters
    ----------
        filepath : str
            Filename or entire path of the eventlist.
            
        frequencies : float or list of floats
            testfrequencies to test with the epochfolding algorithm.
            If a single float, a narrow range around the given frequency will be performed.
            For epochfolding with a single frequency, please use epochfolding_single().
            
    Optional Parameters
    -------------------
        nbin : int
            Number of bins in the folded eventlist.
            Defaults to 32.
            
        oversampling : int
            Oversampling of the folded eventlist.
            Defaults to 10.
            To be specified when `frequencies` is a single float.
            Together with `obs_length` determines the precision of the testfrequencies.
        
        obs_length : float
            Observation length
            To be specified when `frequencies` is a single float.
            Together with `oversampling` determines the precision of the testfrequencies.
        
        plot : bool
            Whether to plot the statistics of the epoch folding search.
            
        true_frequency : float
            True frequency of the signal (if known).
            
    Returns
    -------
        frequencies : array-like
            epoch folding trial frequencies.
            
        effreq : array-like
            frequency grid of the folded timeseries.
            
        efstat : array-like
            epoch folding statistics (chi2-values) for the frequency bins in `effreq`.
    
    """
    events = EventList().read(filepath, format)
        
    if isinstance(frequencies, float): 
        if obs_length is None:
            raise ValueError(
                "Please specify the observation length if you enter only one testfrequency."
            )
            
        df_min = 1/obs_length
        df = df_min / oversampling
        frequencies = np.arange(true_frequency - number_testf * df, true_frequency + (number_testf+1) * df, df)
    elif isinstance(frequencies, Iterable):
        frequencies = frequencies
    else:
        raise ValueError(
            "testfrequencies should be a float or a list of floats!"
        )
    #print(events.time, frequencies, nbin)
    #print(type(events.time[0]), type(frequencies[0]), type(nbin))
    print(frequencies)
    effreq, efstat = epoch_folding_search(events.time, frequencies, nbin=nbin)
    
    if plot:
        plot_efstat(effreq, efstat, nbin, true_frequency=true_frequency)
    
    if save:
         with h5py.File(output, 'w') as output_file:
            savehdf5(effreq, efstat, output_file)
        
    return frequencies, effreq, efstat


def plot_efstat(effreq, efstat, nbin, output='chi2', label='EF statistics', true_frequency=None):
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
            
        true_frequency : float
            True frequency of the signal (if known).
            Default is None.
        
    """
    plt.figure(figsize=figsize)
    plt.plot(effreq, efstat, label=r'$\chi^2$ landscape', color='blue')
    plt.axhline(nbin - 1, ls='--', lw=2, color='k', label=r'd.o.f. $ = N_\mathrm{{bin}} - 1 = {}$'.format(nbin-1))
    plt.axvline(effreq[np.argmax(efstat)], lw=2, alpha=0.5, color='grey', label='$f_\mathrm{{obs}} = {} \,\mathrm{{Hz}}$'.format(round(effreq[np.argmax(efstat)], 3)))
    if isinstance(true_frequency, float):
        plt.axvline(true_frequency, lw=2, alpha=0.5, color='r', label=r'$f_\mathrm{true} = {' + str(true_frequency) + '\,\mathrm{Hz}}$')
    plt.xlabel(r'Frequency $f$ [Hz]')
    plt.ylabel(r'$\chi^2$ test statistic')
    plt.legend()
    plt.savefig(output + '.pdf', bbox_inches='tight')
    

def fit_stats( frequencies, nbin, effreq, efstat, obs_length, pulseshape='arbitrary', plot=False, true_frequency=None ):
    """Fits a gaussian or a sinc-function to the epoch folding statistics around ***the first peak***.
    
    Parameters
    ----------
        frequencies : array-like
            epoch folding trial frequencies.
            
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
            
        true_frequency : float
            True frequency of the signal (if known).
            
    Returns
    -------
        cand_freq_ef : array-like
            array containing the x position of the peaks above threshold. If no
            peaks are above threshold, an empty list is returned. The array is
            sorted by inverse value of statistics.
            
        cand_stat_ef :
            for each peak in `cand_freq_ef`, give the corresponding chi2 value. Empty if no peaks
            above threshold.
            
        bestfit : function
            The best-fit function, accepting x as input
            and returning the best-fit model as output.
        
    """
    
    # search best candidate frequencies
    df_min = 1/obs_length
    ntrial = (frequencies[-1] - frequencies[0]) / df_min
    ef_detlev = fold_detection_level(nbin, epsilon=0.001, ntrial=len(effreq))
    cand_freqs_ef, cand_stat_ef = search_best_peaks(effreq, efstat, ef_detlev)
    
    if pulseshape == 'arbitrary':
        bestfit = fit_gaussian(effreq, efstat-(nbin-1),amplitude=max(efstat-(nbin-1)), 
                mean=cand_freqs_ef[0], stddev=1/(np.pi*obs_length))
    elif pulseshape == 'sine':
        bestfit = fit_sinc(effreq, efstat-(nbin-1),amp=max(efstat-(nbin-1)), mean=cand_freqs_ef[0], 
            obs_length=obs_length)
    else:
        raise ValueError(
            "pulseshape should be either `arbitrary` or `sine`!"
        )
        
    if plot:
        plot_stats_fit(bestfit, effreq, efstat, nbin, output='chi2_fit', label='EF statistics', true_frequency=true_frequency)
    
    return cand_freqs_ef, cand_stat_ef, bestfit

def plot_stats_fit(bestfit, effreq, efstat, nbin, output='chi2_fit', label='EF statistics', true_frequency=None):
    """Plots the chi2 distribution of the epochfolding scan in frequency space.
    
    Parameters
    ----------
        bestfit : function
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
            
        true_frequency : float
            True frequency of the signal (if known).
            Default is None.
        
    """
        
    plt.figure(figsize=(15, 5))
    plt.plot(effreq, efstat-(nbin-1), label=label, color='blue')
    plt.plot(effreq, bestfit(effreq), label='Best fit', color='red', ls='--')
    if isinstance(true_frequency, float):
        plt.axvline(true_frequency, alpha=0.5, color='r', label='Correct frequency $f_\mathrm{true} = {' + str(true_frequency) + '}$')
    plt.axvline(bestfit.mean[0], alpha=0.5, label=r'Fit frequency $f_\mathrm{fit} = $' + '{}'.format(bestfit.mean[0]))
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
    print(arguments)
    filepath = arguments['-i']
    mode = arguments['--mode']
    nbin = int(arguments['--nbin'])
    oversampling = int(arguments['--oversampling'])
    obs_length = float(arguments['--obs_length'])
    true_frequency = float(arguments['--true_frequency'])
    if arguments['--frequencies'] is None:
        #obs_length = 1000
    #    oversampling = 10
    #    df_min = 1/obs_length
    #    df = df_min / oversampling
    #    frequencies = np.arange(true_frequency - 100 * df, true_frequency + 100 * df, df)
    #    print(len(frequencies))
        df_min = 1/obs_length
        df = df_min / oversampling
        frequencies = np.arange(true_frequency - float(arguments['--number_testf']) * df, true_frequency + float(arguments['--number_testf']) * df, df)
    elif arguments['--frequencies'][0] != '[':
        frequencies = float(arguments['--frequencies'])
    elif arguments['--frequencies'][0] == '[':
        #print(arguments['--testperiods'].strip('[]').split(','))
        frequencies = list( map(float, arguments['--frequencies'].strip('[]').split(',')) )
    plot = arguments['--plot']
    save = arguments['--save']
    pulseshape = arguments['--pulseshape']
    
    if mode == 'single':
        if arguments['--frequencies'] is None:
            epochfolding_single(filepath, frequencies, output='profile', plot=plot, save=save, outputdir=arguments['--outputdir'])
        else:
            if nbin is None:
                 epochfolding_single(filepath, frequencies, output='profile', plot=plot, save=save)
            else:
                epochfolding_single(filepath, frequencies, nbin=nbin, output='profile', plot=plot, save=save)
        
    elif mode == 'scan':
        frequencies, effreq, efstat = epochfolding_scan(filepath, true_frequency, nbin=nbin, oversampling=oversampling, 
                          obs_length=obs_length, plot=plot, true_frequency=true_frequency, save=save, number_testf=float(arguments['--number_testf']))
        #print(fit)
        if arguments['--fit']:
            fit_stats( frequencies, nbin, effreq, efstat, obs_length, pulseshape=pulseshape, plot=plot, true_frequency=true_frequency )
            
    elif mode == 'all':
        epochfolding_single(filepath, frequencies, nbin=nbin, output='profile', plot=plot, save=save)
        
        frequencies, effreq, efstat = epochfolding_scan(filepath, true_frequency, nbin=nbin, oversampling=oversampling, obs_length=obs_length, plot=plot, true_frequency=true_frequency)
        
        cand_freqs_ef, cand_stat_ef, fit = fit_stats( frequencies, nbin, effreq, efstat, obs_length,
                                  pulseshape=pulseshape, plot=plot, true_frequency=true_frequency )
        