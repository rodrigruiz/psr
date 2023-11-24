import numpy as np
import awkward as ak

from scipy.stats import chi2



def ak_median(array):
    """Calculates the median on axis=1 of an awkward array.
    Based on: https://github.com/scikit-hep/awkward/issues/741
    
    Parameters
    ----------
    array : array
        Input array.
        
    Returns
    -------
    median : array
        Medians of axis=1.
    """
    
    sorted_array = ak.sort(array, axis=1)
    sorted_without_empties = sorted_array.mask[ak.num(sorted_array) != 0]
    
    
    low_index = ak.values_astype(np.floor((ak.num(sorted_array, axis=1) - 1) / 2), np.int64)
    high_index = ak.values_astype(np.ceil((ak.num(sorted_array, axis=1) - 1) / 2), np.int64)
    
    low = sorted_without_empties[ak.from_regular(low_index[:, np.newaxis])]    
    high = sorted_without_empties[ak.from_regular(high_index[:, np.newaxis])]
    
    median = ((low + high) / 2)[:, 0]
    
    return median



def ak_RMedianSE(array, ddof=0):
    """Calculates the Root Median Square Error.
    As an equivalent to the standard deviation corresponding to the median instead of the mean.
    
    Parameters
    ----------
    array : array
        Input array.
    ddof : int
        Delta Degrees of Freedom. The divisor used in calculations is N - ddof, where N represents the number of elements. By default ddof is 0.
    Returns
    -------
    median : array
        Medians of axis=1.
    """
    median = ak_median(array)
    
    RMedianSD = np.sqrt(ak.sum( (median - array)**2 , axis=1) / (ak.count(array, axis=1)-ddof))
    
    return RMedianSD



def RedNoiseFilter(fft_spectrum, indices_or_sections, ddof=0, central_tendency='mean', exclude_max=False, return_segments=False):
    """Whitens the given frequency spectrum.
    
    Parameters
    ----------
    fft_abs : array
        The spectrum to normalize (without 0'th frequency bin). 
    indices_or_sections : int or array
        If indices_or_sections is an integer, N, the array will be divided into N equal arrays along axis. If such a split is not possible, an error is raised.
        If indices_or_sections is a 1-D array of sorted integers, the entries indicate where along axis the array is split.
    ddof : int, optional
        Means Delta Degrees of Freedom. The divisor used in calculations is N - ddof, where N represents the number of elements. By default ddof is 0.
        #In the current implementation only available for central_tendency='mean'.
    exclude_max : bool, optional
        If True, the maximum value value of each segment is excluded in the calculation of the normalization parameters. Default is False.
    return_segments : bool, optional
        Additionaly returns segments if True for further statistical analysis. Default is False.
        
    Returns
    -------
    normalized_fft_spectrum : array
        An array containing the normalised frequency spectrum.
    normalized_segments : array (Awkward-Array)
        An awkward array containing the normalized segments.
    """    
    
    # Split spectrum into smaller segments
    segments = ak.Array(np.array_split(fft_spectrum, indices_or_sections))
    # Another possible approach: ak.unflatten(); requires number of elements in each sub-array, problems appear, the given number of elements in all sub-arrays was smaller thant the number of array in the input array.
    
    # Exclude the highes value of each segment for the mean/median and std calculation to decrease the influence of a potential signal on the normalization process
    if exclude_max == True:
        idx = segments != ak.max(segments, axis=1)
    elif exclude_max == False:
        idx = segments <= ak.max(segments, axis=1)
    
    if central_tendency=='mean':
        mean = ak.mean(segments[idx], axis=1)
        std = ak.std(segments[idx], axis=1, ddof=ddof)

        normalized_segments = (segments - np.array(mean)) / np.array(std)
        # mean and std are wrapped with np.array() because normalized_segments has then the type '5 * var * float64' and not '5 * option[var * float64]'. I'm not aware if any of these two options causes any problems, but the first one seems nicer.

    if central_tendency=='median':
        median = ak_median(segments[idx])
        std = ak_RMedianSE(segments[idx], ddof=ddof)
        normalized_segments = (segments - np.array(median)) / np.array(std)

    normalized_fft_spectrum = ak.to_numpy(ak.flatten(normalized_segments))

    
    
    if return_segments==False:
        return normalized_fft_spectrum 
    
    elif return_segments==True:
        return normalized_fft_spectrum, normalized_segments


    
def loc_RNF(df):
    """Get the location parameter of the chi-squared distribution created by the RedNoiseFilter.
    
    Parameters
    ----------
    df : int
        Degrees of freedom of the chi-squared distribution.
    
    Returns
    -------
    loc_RNF : float
        Location parameter of the chi-squared distribution.
    
    """
    
    loc_RNF = - np.sqrt(df / 2)
    return loc_RNF


def scale_RNF(df):
    """Get the scale parameter of the chi-squared distribution created by the RedNoiseFilter.
    
    Parameters
    ----------
    df : int
        Degrees of freedom of the chi-squared distribution.
    
    Returns
    -------
    scale_RNF : float
        Scale parameter of the chi-squared distribution.
    
    """
    
    scale_RNF = 1 / np.sqrt(2 *df)
    return scale_RNF
    

def chi2_RNF(df):
    """Get the chi-square distribution created by the RedNoiseFilter.
    
    Parameters
    ----------
    df : int
        Degrees of freedom of the chi-squared distribution.
    
    Returns
    chi2_RNF : scipy.stats._distn_infrastructure.rv_continuous_frozen
        Frozen scipy.stats chi-squared distribution.
    
    """
    
    chi2_RNF = chi2(df, loc=loc_RNF(df), scale=scale_RNF(df))
    return chi2_RNF