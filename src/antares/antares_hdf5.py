"""
Author: Maximilian Eff
Date: 2022-12-16
Description: This is a library to read hdf5 files containing the module rates recorded by the ANTARES telescope.
"""
import numpy as np
import h5py as h5py

from astropy import units as u
from astropy.coordinates import EarthLocation
from astropy.time import Time, TimeDelta


# Collection of functions, to ectract data from ANTARES HDF5 files


def get_start_time(run):
    """Get start time of a run as astropy's Time object.
    
    Parameters
    ----------
    run : HDF5 file
        HDF5 file object of a run.
        
    Returns
    -------
    start_time : astropy Time
        Start time of the run.
    """
    
    # ANTARES times assumed to be unix time. [citation needed]
    # How are leap seconds handled?
    
    # ANTARES Location
    # Reference: https://doi.org/10.1016/j.nima.2011.06.103 - ANTARES: The first undersea neutrino telescope
    # Footnote 6: 42°48N, 6°10E, 2475 m depth
    # Section 2.2.1.: 100 m from the seabed, 25 OMFs with distance of 14.5 m
    # => Center at half of the 24 gaps between the 25 storey
    height = -2475*u.meter + 100*u.meter + (25-1)/2 * 14.5*u.meter
    antares_location = EarthLocation.from_geodetic(lat='42°48N', lon='6°10E', height=height)

    start_time = Time(run['header/gps_start_time_unix'][()].decode(), format='unix', location=antares_location)
    return start_time



def get_timeslice_duration(run):
    """Get timeslice duration of a run as astropy's TimeDelta object.
    
    Parameters
    ----------
    run : HDF5 file
        HDF5 file object of a run.
    
    Returns
    -------
    timeslice_duration : astropy TimeDelta
        Timeslice duration of the run.
    """
     
    timeslice_duration = TimeDelta(np.double(run['header/timeslice_duration_in_millisecs'][()].decode())*1e-3, format='sec')
    return timeslice_duration



def get_num_timeslices(run):
    """Get number of timeslices of a run.
    
    Parameters
    ----------
    run : HDF5 file
        HDF5 file object of a run.
    
    Returns
    -------
    num_timeslices : int
        Numer of timeslices in the run.
    """
    
    #Needs to be first converted to float, before int, as the string end with a '.0'.
    num_timeslices = int(np.double(run['header/nb_of_timeslices_in_run'][()].decode()))
    return num_timeslices



def get_lcmID_dict(run):
    """Get index dictionary of lcmIDs for a run.
    
    Parameters
    ----------
    run : HDF5 file
        HDF5 file object of a run.
    
    Returns
    -------
    lcmID_dict : dict
        Dictionary of lcmIDs (keys) and corresponding ars_rates array indices (values) for the run.
    """
    
    keys = run['lcmID_dict/keys'][()]
    values = run['lcmID_dict/values'][()]
    lcmID_dict = dict(zip(keys, values))
    return lcmID_dict



def get_num_active_ars(run):
    """Get number of active ars for a run.
    
    Parameters
    ----------
    run : HDF5 file
        HDF5 file object of a run.
    
    Returns
    -------
    num_active_ars : int
        Number of active ars of the run.
    """

    num_active_ars = run['num_active_ars'][()]
    return num_active_ars



def get_total_rates(run):
    """Get total rates for a run.
    
    Parameters
    ----------
    run : HDF5 file
        HDF5 file object of a run.
    
    Returns
    -------
    total_rates : array
        Total detector rates time series of rateOff and rateOn.
    """

    if 'total_rates' in run.keys():
        total_rates = run['total_rates'][()]
    else:
        total_rates = run['rates/total_rates'][()]
    return total_rates



def get_ars_rates(run, source_sel=None):
    """Get ars rates of run.
    
    Parameters
    ----------
    run : HDF5 file
        HDF5 file object of a run.
    source_sel : slice , optional
        Slices of the HDF5 ars array that should be copied (default is None). (Use np.s_)
    
    Returns
    -------
    ars_rates : array
        Selected ars rates of the run.
    """

    # Can be probably further improved. Giving no input for source_sel apparently takes slightly longer that giving np.s_[:,:,:,:]
    # The feature source_sel, that would allow to only extract a subarray is not necessary and could also be abandond, if not required
    # Constructing ars_rates array two times is maybe unnecessary
    if 'rates' not in run.keys():
        raise Exception('File contains no ars rates.')        
    ars_rates = np.zeros(run['rates/ars_rates'].shape)   
    if source_sel != None:
        ars_rates = np.zeros(ars_rates[source_sel].shape, order='C')     
    run['rates/ars_rates'].read_direct(ars_rates, source_sel=source_sel, dest_sel=None)
    return ars_rates




