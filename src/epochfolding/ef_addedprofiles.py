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

#from generate_split_timeseries import generate_pulse_train_gauss, gauss
#from stingray import Lightcurve
from stingray.pulse.pulsar import profile_stat #for old stingray version
#from stingray.pulse.pulsar import ef_profile_stat for new stingray version, change function name in line117 or import as profile_stat
#"""
from stingray.events import EventList
from stingray.pulse.search import _folding_search, _profile_fast
from stingray.utils import jit, HAS_NUMBA
#from stingray.pulse.pulsar import ef_profile_stat
from stingray.stats import fold_detection_level
from stingray.pulse.search import plot_profile, search_best_peaks
from stingray.pulse.modeling import fit_gaussian, fit_sinc
#"""
from epochfolding.add_folded_profiles import add_profiles#, add_folded_profiles
from epochfolding.stingray_epochfolding import plot_efstat

def save_to_ascii(frequencies, chi2, output_dir='./', output='chi2_added_'):
    """Helping function.
       Saves the output of fold_events() to an ascii file.
    """

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        
    data = Table()
    data['frequency'] = frequencies
    data['chi2'] = chi2
    
    ascii.write(data, output_dir + output + '.dat', format='ecsv', overwrite=True)

def epoch_folding_search(times, frequencies, input_dir, filepattern, outputdir, nbin=128, segment_size=np.inf, expocorr=False, gti=None, weights=1, fdots=0):
    """Performs epoch folding at trial frequencies in photon data.

    If no exposure correction is needed and numba is installed, it uses a fast
    algorithm to perform the folding. Otherwise, it runs a *much* slower
    algorithm, which however yields a more precise result.
    The search can be done in segments and the results averaged. Use
    segment_size to control this

    Parameters
    ----------
    times : array-like
        the event arrival times

    frequencies : array-like
        the trial values for the frequencies

    Other Parameters
    ----------------
    nbin : int
        the number of bins of the folded profiles

    segment_size : float
        the length of the segments to be averaged in the periodogram

    fdots : array-like
        trial values of the first frequency derivative (optional)

    expocorr : bool
        correct for the exposure (Use it if the frequency is comparable to the
        length of the good time intervals). If True, GTIs have to be specified
        via the ``gti`` keyword

    gti : [[gti0_0, gti0_1], [gti1_0, gti1_1], ...]
        Good time intervals

    weights : array-like
        weight for each time. This might be, for example, the number of counts
        if the times array contains the time bins of a light curve

    Returns
    -------
    (fgrid, stats) or (fgrid, fdgrid, stats), as follows:

    fgrid : array-like
        frequency grid of the epoch folding periodogram
    fdgrid : array-like
        frequency derivative grid. Only returned if fdots is an array.
    stats : array-like
        the epoch folding statistics corresponding to each frequency bin.
    """
    
    if expocorr or not HAS_NUMBA or isinstance(weights, Iterable):
        if expocorr and gti is None:
            raise ValueError("To calculate exposure correction, you need to" " specify the GTIs")

        def stat_fun(t, f, fd=0, **kwargs):
            return profile_stat(add_profiles(input_dir, filepattern, f, outputdir)[1])
        
        return _folding_search(
            stat_fun,
            times,
            frequencies,
            segment_size=segment_size,
            use_times=True,
            expocorr=expocorr,
            weights=weights,
            gti=gti,
            nbin=nbin,
            fdots=fdots,
        )
    
    def stat_fun(t, f, fd=0, **kwargs):
        return profile_stat(add_profiles(input_dir, filepattern, f, outputdir)[1])
    
    return _folding_search(
            stat_fun,
            times,
            frequencies,
            segment_size=segment_size,
            use_times=True,
            expocorr=expocorr,
            weights=weights,
            gti=gti,
            nbin=nbin,
            fdots=fdots,
        )