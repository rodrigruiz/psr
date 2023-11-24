# This are versions of the barycentric_correction and fill_TimeSeries_gaps functions, that were created for astropy.timeseries.TimeSeries instead of astropy.timeseries.BinnedTimeSeries. They are not needed anymore, but keept here for reference.



def barycentric_correction(timeseries, timeslice_duration, skycoord):
    """Get the brycentric corrected timeseries.
    
    Parameters
    ----------
    timeseries : astropy.timeseries.TimeSeries
        TimeSeries in the earth frame of reference
    timeslice_duration : astropy.time.TimeDelta
        Time difference between two consecutive sample points
    
    Returns
    -------
    TS_bar_cor : astropy.timeseries.TimeSeries
        TimeSeries in the barycentric frame of reference
    
    """
    
    ##### Location of ANTARES currently needs to be manually added to the TimeSeries as the Antares_HDF5.get_start_time() function has the ANTARES location yet not implemented.
    antares_location = EarthLocation.from_geodetic(lat=42.798916, lon=6.1657)
    
    # Calculate the light travel time correction
    dt = timeseries['time'].light_travel_time(skycoord)#, location=antares_location)

    # Add calculated time correction to timeseries
    TS_bar_cor = TimeSeries(time=timeseries['time'] + dt, data={key: timeseries[key] for key in timeseries.keys() if key != 'time'})
     
    return TS_bar_cor






def fill_TimeSeries_gaps(timeseries, timeslice_duration):
    """Get a continuous TimeSeries with equidistant samples. Gaps are filled with the corresponding mean.
    
    Parameters
    ----------
    timeseries : astropy.timeseries.TimeSeries
        TimeSeries with gaps or samples to close together.
    timeslice_duration : astropy.time.TimeDelta
        Intendet time difference between two consecutive sample points
    
    Returns
    -------
    TS_filles : astropy.timeseries.TimeSeries
        Evenly time spaced TimeSeries
    """
    
    # Get start end end of timeseries
    ## Using min/max instead of timeseries['time'][0/-1] has the advantage, that timeseries does not need to be sorted
    start_time = timeseries['time'].min()
    end_time = timeseries['time'].max()


    # Calculate the number of samples, the new TimeSeries requires
    n_samples = int(np.round( (end_time + timeslice_duration - start_time) / timeslice_duration ))
    
    # Create a new empy TimeSeries
    TS_filled = TimeSeries(time_start=start_time, time_delta=timeslice_duration, n_samples=n_samples, data={key: np.zeros(n_samples) for key in timeseries.keys() if key != 'time'})
    
    # Calculate indices into wich the values of timeseries will be filled
    ## np.round() puts the input bin to its closest continuous bin
    ## As the input is usually binned data, maybe it would be better to insead use np.ceil or np.floor, as this would always either round up or down.
    idx = np.array(np.round( (timeseries['time'] - start_time) / timeslice_duration ), dtype=int)
    
    # Calculate mask that covers all entries in idx
    mask = np.ones(n_samples, dtype=bool)
    mask[idx] = False
    
    # Fill new TimeSeries with values from timeseries
    for key in timeseries.keys():
        if key == 'time':
            continue
        else:
            # Fill values of timeseries into new TimeSeries
            # If one index is containted more than once in idx (possible due to barrycentric correction), the last one will be used, as the previous ones get overwritten
            TS_filled[key][idx] = timeseries[key]
            # Fill gaps with the mean of the timeseries column
            TS_filled[key][mask] = np.mean(timeseries[key])
    
    return TS_filled