import numpy as np
from astropy.time import Time, TimeDelta
from astropy.timeseries import BinnedTimeSeries, TimeSeries
from astropy.coordinates import EarthLocation
import astropy.units as u
import h5py
from plens.TimeSeries import antares_location
from plens.PulseModel import MVMD, sinusoid
from stingray import EventList, Lightcurve

def readEventList(file):
    """Read EventList from HDF5-file constucted by CreateEventlist. 
    
    Parameters
    ----------
    file : h5py.File
        Input file to read EventList.
    
    Returns
    -------
    timeseries : astropy.timeseries.BinnedTimeSeries
        The BinnedTimeSeries has four columns: 'time_bin_start'
    
    """
    #with h5py.File(file) as f:
    #print(file.keys())
   # timeseries = TimeSeries(time=Time(file['timeseries/time'][()], format='unix'))
    timeseries = TimeSeries.read(file, format='hdf5', time_column='time', time_format='unix')
    """
    timeseries = BinnedTimeSeries(time_bin_start=Time(file['timeseries/time_bin_start'][()], 
                                                      format='unix', location=antares_location()), 
                                  time_bin_size=TimeDelta(file['timeseries/time_bin_size'][()], format='sec'), 
                                  data={'rateOff': file['timeseries/rateOff'][()]*1e3*u.Hz, 
                                        'rateOn':  file['timeseries/rateOn' ][()]*1e3*u.Hz
                                        }
                                  )
    """
    #timeslice_duration = TimeDelta(file['timeseries/timeslice_duration'], format='sec')
    
    return timeseries

def saveEventList(timeseries, file):
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
    #print(timeseries['time'])
    file['timeseries/time'] = timeseries['time'].to_value('unix')
    #file['timeseries/time_bin_size'] = timeslice_duration.to_value('sec')
    #file['timeseries/rateOff'] = timeseries['rateOff'].to_value()
    #file['timeseries/rateOn'] = timeseries['rateOn'].to_value()
    #file['timeseries/timeslice_duration'] = timeslice_duration.to_value('sec')
    # time_bin_size is for the variable-size bin size, due to barycentric correction
    # timeslice_duration is the underlying constant ANTARES timeslice duration
    
    return
    

def barycentric_correction(timeseries, skycoord):
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
    dt = timeseries.time.light_travel_time(skycoord=skycoord, kind='barycentric', location=antares_location())
    # Calculate the light travel time correction for the end of the last bin
    #dt_lb = timeseries.time_bin_end[-1].light_travel_time(skycoord=skycoord, kind='barycentric', location=antares_location())
    
    # Calculate new time_bin_start
    time_bin_start = timeseries.time.tdb + dt
    
    # Calculate new time_bin_end
    # Due to the correction, the time_bin_size is not constant anymore
    # Use the start of the following bin as end time
    #time_bin_end = Time([time_bin_start[1:], timeseries.time_bin_end[-1].tdb + dt_lb])

    # Add calculated time correction to timeseries
    
    TS_bar_cor = TimeSeries(time=time_bin_start)
                                  #time_bin_end=time_bin_end,
                                  #data={key: timeseries[key] for key in timeseries.keys() 
                                  #      if key not in ['time_bin_start', 'time_bin_size']})
    
    return TS_bar_cor


def injectSignal( time, bin_time, pulseshape, frequency, baseline, a, phi, kappa=None ):
    #
    # Hannes: Redistribute events rather than just add new ones at some point!
    #
    """Injects a signal with a MVM pulseshape into an existing time sequence (eventlist).
    
    Parameters
    ----------
        time : np.array
        
        bin_time : float
        
        frequency : float
            Frequency of the pulse train.
            
        baseline : float
            Offset along the y-axis.
            
        a : float
            Amplitude of the Pulse. Equates to the area of one pulse.
            
        phi : float
            Phase offset of the pulse train.
            
        kappa : float
            Shape parameter giving the width of the function.
        
    Returns
    -------
        np.array
        New times with injected pulsetrain.
        
    """
    
    if pulseshape == 'mvm':
        counts = MVMD(time, frequency, phi, kappa, a, baseline=baseline)
    elif pulseshape == 'sine':
        counts = sinusoid(time, frequency, baseline, a, phi)

    lc = Lightcurve(time, counts, dt=bin_time, skip_checks=True)

    ev = EventList()
    ev.simulate_times(lc)

    return np.sort(np.concatenate((time, ev.time)))


def injectSignalRedistribute( time, bin_time, pulseshape, frequency, baseline, a, phi, kappa=None ):

    """Injects a signal with a MVM pulseshape into an existing time sequence (eventlist) by first deleting random counts to ensure an unchanged total count number.
    
    Parameters
    ----------
        time : np.array
        
        bin_time : float
        
        frequency : float
            Frequency of the pulse train.
            
        baseline : float
            Offset along the y-axis.
            
        a : float
            Amplitude of the Pulse. Equates to the area of one pulse.
            
        phi : float
            Phase offset of the pulse train.
            
        kappa : float
            Shape parameter giving the width of the function.
        
    Returns
    -------
        np.array
        New times with injected pulsetrain.
        
    """
    
    if pulseshape == 'mvm':
        counts = MVMD(time, frequency, phi, kappa, a, baseline=baseline)
    elif pulseshape == 'sine':
        counts = sinusoid(time, frequency, baseline, a, phi)

    # Create a light curve with the desired signal
    lc = Lightcurve(time, counts, dt=bin_time, skip_checks=True)
    
    # Simulate event times from the light curve
    ev = EventList()
    ev.simulate_times(lc)
    
    new_event_times = ev.time
    
    # Determine the number of new events
    num_new_events = len(new_event_times)
    
    # Randomly select events to remove from the original times
    if num_new_events > 0:
        indices_to_remove = np.random.choice(len(time), num_new_events, replace=False)
        remaining_original_events = np.delete(time, indices_to_remove)
    else:
        remaining_original_events = time
    
    # Combine the remaining original events with the new events
    combined_events = np.sort(np.concatenate((remaining_original_events, new_event_times)))

    return combined_events