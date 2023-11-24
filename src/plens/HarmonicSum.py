import numpy as np
from scipy.stats import chi2
from numba import jit

def HarmonicSum(spectrum, H):
    """Harmonic Sum
        Calculates the harmonic sum of H harmonics.
    
    Parameters
    ----------
        spectrum : array
            Spectrum (without 0th frequency bin) to be harmonic summed.
        H : int
            Total number of harmonics summed.
    
    Returns
    -------
        HS_spectrum : array
            Harmonic summed spectrum. Returns only the valid region for H.
            len(HS_spectrum) == len(spectrum)//H
    """
    
    # Create a dummy spectrum with 0th frequency bin, this simplifies the index slicing later
    dummy_spectrum = np.insert(spectrum, 0, 0)
    
    # Create Harmonic Sum output array, length of the valid region for harmonic summing of order H
    HS_spectrum = np.zeros(len(dummy_spectrum[1:])//H)
    
    # Accelerate this loop using numba
    @jit(nopython=True)
    def HSum(HS_spectrum, dummy_spectrum):
        # Loop over each to-be filled element in HS
        for h in range(1,len(HS_spectrum)+1):
            HS_spectrum[h-1] = np.sum(dummy_spectrum[h:h*H+1:h])
        return
    
    HSum(HS_spectrum, dummy_spectrum)

    return HS_spectrum

def loc_HS(df, H):
    """Get the location parameter of the chi-squared distribution created by the RedNoiseFilter.
    
    Parameters
    ----------
    df : int
        Degrees of freedom of the chi-squared distribution.
    H : int
        Highes harmonic order.
    
    Returns
    -------
    loc_HS : float
        Location parameter of the chi-squared distribution.
    
    """
    
    loc_HS = - H * np.sqrt(df / 2)
    return loc_HS


def scale_HS(df):
    """Get the scale parameter of the chi-squared distribution created by the Harmonic Summing.
    
    Parameters
    ----------
    df : int
        Degrees of freedom of the chi-squared distribution.
    
    Returns
    -------
    scale_HS : float
        Scale parameter of the chi-squared distribution.
    
    """
    
    scale_HS = 1 / np.sqrt(2 * df)
    return scale_HS
    

def chi2_HS(df, H):
    """Get the chi-square distribution created by the Harmonic Summing.
    
    Parameters
    ----------
    df : int
        Degrees of freedom of the chi-squared distribution.
    H : int
        Highes harmonic order.
    
    Returns
    chi2_HS : scipy.stats._distn_infrastructure.rv_continuous_frozen
        Frozen scipy.stats chi-squared distribution.
    
    """
    
    chi2_HS = chi2(df*H, loc=loc_HS(df,H), scale=scale_HS(df))
    return chi2_HS