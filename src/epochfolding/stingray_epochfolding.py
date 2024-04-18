import os, glob
import h5py
from docopt import docopt
from astropy.table import Table
from astropy.io import ascii
from collections.abc import Iterable
import numpy as np
import matplotlib.pyplot as plt
#plt.style.use('~/software/Psr/src/latex.mplstyle')
import seaborn as sb
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (10, 6)
#import astropy.time.core.Time
from astropy.time import Time
#from generate_split_timeseries import generate_pulse_train_gauss, gauss
#from stingray import Lightcurve
from stingray.events import EventList
from stingray.pulse.search import epoch_folding_search, _folding_search
from stingray.pulse.pulsar import fold_events
from stingray.stats import fold_detection_level
from stingray.pulse.search import plot_profile, search_best_peaks
from stingray.pulse.modeling import fit_gaussian, fit_sinc

def load_eventlist(directory):
    file_pattern = os.path.join(directory, 'eventlist_*.dat')
    eventlist_files = glob.glob(file_pattern)
    
    if not eventlist_files:
        print(f"No timeseries files found in the directory: {directory}")
        return
    
    key = lambda x: int(x.split('_')[-1].split('.')[0])
    eventlist_files.sort(key=key)
    
def epochfolding_single(filepath, frequencies, nbin=32, expocorr=False, gti=None, output='profile', plot=True, save=True, format='ascii', outputdir='./folded_events/'):
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
    
    events = EventList().read(filepath, format)
    #print(type(events.time), type(events.time[0]))
    _ef_single(events, frequencies, nbin=nbin, expocorr=expocorr, gti=gti, output=output, plot=plot, save=save, outputdir=outputdir, root_dir=None, format=format)

def _ef_single(events, frequencies, nbin=32, expocorr=False, gti=None, output='profile', plot=True, save=True, outputdir='./folded_events/', root_dir=None, format='ascii'):
    """Helping function.
       Calls _fold_events and plots the folded events (if plot).
    """
    
    if plot:
        fig, ax = plt.subplots()
    
    # perform epochfolding and plot resulting profile
    if isinstance(frequencies, Iterable):
        
        label = []
        for p in frequencies:
            phase_bins, profile, profile_err = _fold_events(events.time, p, nbin=nbin, expocorr=expocorr, gti=gti, save=save, output=output+f'_{p}', outputdir=outputdir, root_dir=root_dir, format=format)
            if plot:
                ax = plot_profile(phase_bins, profile, ax=ax)
        
                label.append(r'$f_\mathrm{test} = ' +  str(p) + '\, \mathrm{Hz}$ ')
                ax.legend(label)

                fig.savefig(output + '.png', bbox_inches='tight')

    elif isinstance(frequencies, float):
        phase_bins, profile, profile_err = _fold_events(events.time, frequencies, nbin=nbin, expocorr=expocorr, gti=gti, save=save, output=output+f'_{frequencies}', outputdir=outputdir, root_dir=root_dir, format=format)
        if plot:
            ax = plot_profile(phase_bins, profile, ax=ax)
            fig.savefig(output + '.png', bbox_inches='tight')
        
    else:
        raise ValueError(
            "frequencies should be either a float of a list of floats!"
        )
    
    plt.close()

def _fold_events(times, frequency, nbin=32, expocorr=False, gti=None, save=False, output='profile', outputdir='./folded_events/', root_dir=None, format='ascii'):
    """Helping function.
       Calculate the folded pulse profile.
    """
    if isinstance(times[0], Time):
        times = np.array([i.value for i in times])

    phase_bins, profile, profile_err = fold_events(times, frequency, nbin=nbin, expocorr=expocorr, gti=gti)
    if save:
        if format=='ascii':
            save_to_ascii(phase_bins, profile, profile_err, frequency, outputdir=outputdir, root_dir=root_dir, output=output)
        if format=='hdf5':
            data = { 'ef_profile/phase_bins' : phase_bins, 'ef_profile/profile' : profile, 'ef_profile/profile_err' : profile_err }
            save_hdf5(outputdir+output+'.hdf5', **data )
        
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


def save_hdf5(filename, **data):
    with h5py.File(filename, 'w') as file:
        for key, value in data.items():
            file[key] = value
    
    return
    
def epochfolding_scan(filepath, frequencies, nbin=32, oversampling=10, number_testf=200, obs_length=None, expocorr=False, gti=None, plot=True, save=False, true_frequency=None, format='ascii', output='chi2s.hdf5'):
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
        
        expocorr : bool
            correct for the exposure (Use it if the period is comparable to the
            length of the good time intervals). If True, GTIs have to be specified
            via the ``gti`` keyword
        
        gti : list of tuples
            Good Time Intervals
            
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
        frequencies = np.arange(frequencies - number_testf/2 * df, frequencies + number_testf/2 * df, df)
    elif isinstance(frequencies, Iterable):
        frequencies = frequencies
    else:
        raise ValueError(
            "testfrequencies should be a float or a list of floats!"
        )
    
    #print(type(events.time), type(events.time[0]))
    effreq, efstat = epoch_folding_search(events.time, frequencies, nbin=nbin, expocorr=expocorr, gti=gti)
    
    if plot:
        plot_efstat(effreq, efstat, nbin, output='chi2_{}'.format(filepath), true_frequency=true_frequency)
    
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
    
    plt.figure()
    plt.plot(effreq, efstat, label=label, color='blue')
    plt.axhline(nbin - 1, ls='--', lw=2, color='k', label=r'd.o.f. $ = N_\mathrm{bin} - 1$')
    plt.axvline(effreq[np.argmax(efstat)], lw=3, alpha=0.5, color='k', label='Observed frequency')
    if isinstance(true_frequency, float):
        plt.axvline(true_frequency, lw=3, alpha=0.5, color='r', label=r'Correct frequency $f_\mathrm{true} = {' + str(true_frequency) + '}$')
    plt.xlabel(r'Frequency $f$ [Hz]')
    plt.ylabel(r'EF Statistics ($\chi^2$)')
    plt.legend()
    plt.savefig(output + '.png', bbox_inches='tight')
    

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
        plt.axvline(true_frequency, alpha=0.5, color='r', label=r'Correct frequency $f_\mathrm{true} = {' + str(true_frequency) + '}$')
    plt.axvline(bestfit.mean[0], alpha=0.5, label=r'Fit frequency $f_\mathrm{fit} = $' + '{}'.format(bestfit.mean[0]))
    #plt.xlim( [fg.mean[0] * 0.95, fg.mean[0] * 1.05] )
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('EF Statistics')
    plt.legend()
    plt.savefig(output + '.png', bbox_inches='tight')
    

def get_testfrequencies(principal_f, number_of_testf, df):
    """Defines the testfrequencies around the principal frequency in the interval 
    [principal_f - number_of_testf/2 * df, principal_f  + number_of_testf/2 * df]
    with a df spacing between the testfrequencies.
    
    Parameters
    ----------
        principal_f : float
            Principal frequency around which testfrequencies will be defined.
        
        number_of_testf : int
            Number of testfrequencies to define.
            
        df : float
            Spacing between testfrequencies.
            
    Optional Parameters
    -------------------
        
        frequencies : list
            List of testfrequencies.
        
    """
    
    frequencies = np.arange(principal_f - number_of_testf/2 * df, principal_f  + number_of_testf/2 * df, df)
    
    return frequencies
