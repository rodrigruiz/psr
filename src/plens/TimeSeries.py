import numpy as np
from astropy.time import Time, TimeDelta
from astropy.timeseries import BinnedTimeSeries
from astropy.coordinates import EarthLocation
import astropy.units as u
import h5py

import plens.antares_hdf5 as antares_hdf5



# Collection of function to create astropy BinnedTimeSeries from ANTARES HDF5 files

# Use get_TimeSeries to construct a TimeSeries for each ANTARES run, combine them, using astropy.vstack. Save the combined TimeSeries is a new HDF5-file.
# Use read_TimeSeries to construct a TimeSeries, for the already combined TimeSeries of several runs.
# Use first barycentric_correction, and afterwards fill_TimeSeries_gaps.


def antares_location():
    """Get ANTARES Location.
    
    Returns
    -------
    antares_location : astropy.coordinates.EarthLocation
        ANTARES Location
    
    """
    
    # ANTARES Location
    # Reference: https://doi.org/10.1016/j.nima.2011.06.103 - ANTARES: The first undersea neutrino telescope
    # Footnote 6: 42°48N, 6°10E, 2475 m depth
    # Section 2.2.1.: 100 m from the seabed, 25 OMFs with distance of 14.5 m
    # => Center at half of the 24 gaps between the 25 storey
    height = -2475*u.meter + 100*u.meter + (25-1)/2 * 14.5*u.meter
    antares_location = EarthLocation.from_geodetic(lat='42°48N', lon='6°10E', height=height)
    
    return antares_location


def get_TimeSeries(file):
    """Get TimeSeries for input run HDF5-file.
    
    Parameters
    ----------
    file : h5py.File
        Input file to create a BinnedTimeSeries from
    Returns
    -------
    timeseries : astropy.timeseries.BinnedTimeSeries
        The BinnedTimeSeries has four collumns: 'time_bin_start', 'time_bin_size', 'rateOff', 'rateOn'
    
    """
    
    # Get normalized rates per OpticalModule by dividing by the number of active ars and multiplying by 2 ars per OM
    rates = antares_hdf5.get_total_rates(file) / antares_hdf5.get_num_active_ars(file) * 2
    
    # Check num of timeslices
    n_bins = len(antares_hdf5.get_total_rates(file))
    if n_bins != antares_hdf5.get_num_timeslices(file):
        print('Warning: num of timeslices "nb_of_timeslices_in_run" wrong in HDF5-file: ' + str(file))
    
    # Create BinnedTimeSeries
    timeseries = BinnedTimeSeries(time_bin_start=antares_hdf5.get_start_time(file),
                                  time_bin_size=antares_hdf5.get_timeslice_duration(file),
                                  n_bins=n_bins,
                                  data={'rateOff': rates[:,0]*1e3*u.Hz, 'rateOn': rates[:,1]*1e3*u.Hz})
    
    
    ########## BUG-Fix ##########
    timeseries['rateOn'] = timeseries['rateOn']/2
    ########## Remove if the hdf5 extraction script calculates the mean of two ars for both rateOff and rateOn ##########
    return timeseries


def saveTimeSeries(timeseries, timeslice_duration, file):
    """Saves TimeSeries into HDF5-File.
    
    Parameters
    ----------
        timeseries : astropy.timeseries.BinnedTimeSeries
            Timeseries to store. Usually after barycentric correction and filled gaps.
        timeslice_duration : astropy.time.TimeDelta
            Time difference between two consecutive sample points (bevore correction).
        file : h5py.File
            Output file, the timeseries is saved to.
    
    """  
    
    file['timeseries/time_bin_start'] = timeseries['time_bin_start'].to_value('unix')
    file['timeseries/time_bin_size'] = timeslice_duration.to_value('sec')
    file['timeseries/rateOff'] = timeseries['rateOff'].to_value()
    file['timeseries/rateOn'] = timeseries['rateOn'].to_value()
    file['timeseries/timeslice_duration'] = timeslice_duration.to_value('sec')
    # time_bin_size is for the variable-size bin size, due to barycentric correction
    # timeslice_duration is the underlying constant ANTARES timeslice duration
    
    return
    

def readTimeSeries(file):
    """Read TimeSeries from HDF5-file constructed by CreatTimeSeries.py. 
    
    Parameters
    ----------
    file : h5py.File
        Input file to read TimeSeries.
    
    Returns
    -------
    timeseries : astropy.timeseries.BinnedTimeSeries
        The BinnedTimeSeries has four collumns: 'time_bin_start', 'time_bin_size', 'rateOff', 'rateOn'
    
    """
    
    timeseries = BinnedTimeSeries(time_bin_start=Time(file['timeseries/time_bin_start'][()], 
                                                      format='unix', location=antares_location()), 
                                  time_bin_size=TimeDelta(file['timeseries/time_bin_size'][()], format='sec'), 
                                  data={'rateOff': file['timeseries/rateOff'][()]*1e3*u.Hz, 
                                        'rateOn':  file['timeseries/rateOn' ][()]*1e3*u.Hz
                                        }
                                  )
    
    timeslice_duration = TimeDelta(file['timeseries/timeslice_duration'], format='sec')
    
    return timeseries, timeslice_duration

##### Don't use. Very memory intensive for unclear reasons
def barycentric_correction(timeseries, timeslice_duration, skycoord):
    """Get the brycentric corrected timeseries.
    
    Parameters
    ----------
    timeseries : astropy.timeseries.BinnedTimeSeries
        BinnedTimeSeries in the earth frame of reference
    timeslice_duration : astropy.time.TimeDelta
        Time difference between two consecutive sample points
    
    Returns
    -------
    TS_bar_cor : astropy.timeseries.BinnedTimeSeries
        BinnedTimeSeries in the barycentric frame of reference
    
    """
        
    # Calculate the light travel time correction
    dt = timeseries.time_bin_start.light_travel_time(skycoord=skycoord, kind='barycentric', location=antares_location())
    # Calculate the light travel time correction for the end of the last bin
    dt_lb = timeseries.time_bin_end[-1].light_travel_time(skycoord=skycoord, kind='barycentric', location=antares_location())
    
    # Calculate new time_bin_start
    time_bin_start = timeseries.time_bin_start.tdb + dt
    
    # Calculate new time_bin_end
    # Due to the correction, the time_bin_size is not constant anymore
    # Use the start of the following bin as end time
    time_bin_end = Time([time_bin_start[1:], timeseries.time_bin_end[-1].tdb + dt_lb])

    # Add calculated time correction to timeseries
    TS_bar_cor = BinnedTimeSeries(time_bin_start=time_bin_start,
                                  time_bin_end=time_bin_end,
                                  data={key: timeseries[key] for key in timeseries.keys() 
                                        if key not in ['time_bin_start', 'time_bin_size']})
    

    
    
    
    return TS_bar_cor



def fill_TimeSeries_gaps(timeseries, timeslice_duration, replace_zeros=True):
    """Get a continuous TimeSeries with equidistant samples. Gaps are filled with the corresponding mean.
    
    Parameters
    ----------
    timeseries : astropy.timeseries.BinnedTimeSeries
        BinnedTimeSeries with gaps or samples to close together.
    timeslice_duration : astropy.time.TimeDelta
        Intendet time difference between two consecutive sample points
    replace_zeros : bool, default True
        If True, also fills zero values of timeseries with the mean.
        If False, only fills gaps in the timeseries with the mean.
    
    Returns
    -------
    TS_filles : astropy.timeseries.BinnedTimeSeries
        Evenly time spaced BinnedTimeSeries
    """
    
    # Get start end end of timeseries
    ## Using min/max instead of timeseries['time'][0/-1] has the advantage, that timeseries does not need to be sorted
    start_time = timeseries['time_bin_start'].min()
    end_time = timeseries['time_bin_start'].max()

    # Calculate the number of samples, the new TimeSeries requires
    n_bins = int(np.ceil( (end_time + timeslice_duration - start_time) / timeslice_duration ))

    # Create a new empy BinnedTimeSeries
    TS_filled = BinnedTimeSeries(
        time_bin_start=start_time, 
        time_bin_size=timeslice_duration, 
        n_bins=n_bins,
        data={key: np.zeros(n_bins, dtype=timeseries[key].dtype) for key in timeseries.keys() 
              if key not in ['time_bin_start', 'time_bin_size']},
        units={key: timeseries[key].unit for key in timeseries.keys() 
               if key not in ['time_bin_start', 'time_bin_size']})
    
    # Calculate indices into wich the values of timeseries will be filled
    ## np.round() puts the input bin to its closest continuous bin start time
    ## Maybe it would be better to instead use np.ceil or np.floor, as this would always either round up or down.
    idx = np.array(np.round( (timeseries['time_bin_start'] - start_time) / timeslice_duration ), dtype=int)
    
    # Calculate mask that covers all entries in idx
    mask = np.ones(n_bins, dtype=bool)
    mask[idx] = False
    
    # Fill new TimeSeries with values from timeseries
    for key in timeseries.keys():
        if key in ['time_bin_start', 'time_bin_size']:
            continue
        else:
            # Fill values of timeseries into new TimeSeries
            # If one index is containted more than once in idx (possible due to barrycentric correction),
            # the last one will be used, as the previous ones get overwritten
            TS_filled[key][idx] = timeseries[key]
            # Fill gaps with the mean of the timeseries column
            TS_filled[key][mask] = np.mean(timeseries[key])
            if replace_zeros:
                TS_filled[key][TS_filled[key] == 0] = np.mean(timeseries[key])
      
    return TS_filled

# Notes to fill_TimeSeries_gaps:
    # As the ANTARES summary slices are binned rates and not sampled rates, a different approach to this resampling process would actually be preferable. Instead of filling the input bins into their closest output bin, it would be better to have a weighted resampling process, such that all input bins overlaping with one output bin are averaged, weighted by their contributing time.
    # Such a weighted resampling is probably necessary for correction methods that correct constant frequency derivatives, as there the timeseries is resampled with unevenly spaced output bins, simmilar to the barycentric correction. As in this frequency derivatives correction the resampling is much stronger, the weighted resampling probably becomes necessary.
    # I got stuck creating such a function already in the theoretical preparations, as one needs to consider that one output bin can overlap with several input bins, but also one output bin can also be just one fraction of one single input bin. I have yet not found a way to do this in a nice way.

