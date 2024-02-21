#!/usr/bin/env python3
"""
Analyze time series data and create a combined plot of folded time series.

Usage:
  analyze_timeseries.py --directory=<directory> --period=<period> --reference_mjd=<reference_mjd>
  analyze_timeseries.py (-h | --help)

Options:
  -h --help                   Show this help message and exit.
  --directory=<directory>     Directory where timeseries files are stored.
  --period=<period>           Period for folding.
  --reference_mjd=<reference_mjd> Reference MJD for folding.
"""

import os
import glob
import numpy as np
from astropy.timeseries import TimeSeries
import matplotlib.pyplot as plt
from docopt import docopt
from astropy import units as u  # Added for Quantity
from astropy.time import Time, TimeDelta  # Added for Time
from astropy.timeseries import aggregate_downsample
from astropy.table import vstack
from astropy.stats import sigma_clipped_stats
from decimal import Decimal

def combine_timeseries(directory, flux_column_name='counts'):
    '''combines the timeseries from different files into a single timeseries.'''
    
    file_pattern = os.path.join(directory, 'timeseries_*.dat')
    timeseries_files = glob.glob(file_pattern)

    if not timeseries_files:
        print(f"No timeseries files found in the directory: {directory}")
        return

    # Sort files by the integer value i in the filename
    timeseries_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
    
    combined_ts = None
    
    for file_path in timeseries_files:
        ts = TimeSeries.read(file_path, format='ascii', time_column='time', time_format='mjd')
        #print( combined_ts )
        if combined_ts is None:
            combined_ts = ts
        else:
            combined_ts = vstack([combined_ts, ts])
    combined_ts.write('combined_timeseries.dat', format='ascii.ecsv', overwrite=True)  
    
    return combined_ts

def plot_chi2_landscape(input_file=None, chi2s=None, periods=None, output='chi2_landscape'):
    '''plots the chi2 landscape over the test_period.'''
    if input_file != None:
        #data = np.genfromtxt(input_file, dtype=float, skip_header=1)
        data = np.loadtxt(input_file, dtype=None, skiprows=1, unpack=True)
        print(data)
        chi2s, periods = data[1], data[0]
        
    plt.figure(figsize=(7, 7))
    for chi2, period in zip(chi2s, periods):
        plt.vlines(x=period, ymin=0., ymax=chi2, color='purple')
    #plt.plot(periods, chi2s, 'o', label=r'$\chi^2$ lanscape', color='purple' )
    plt.xlabel('Test periods [days]')
    plt.ylabel(r'$\chi^2$')
    plt.title(r'$\chi^2$ lanscape')
    plt.grid(True)
    #plt.legend()
    plt.savefig( output + '.png' )
    plt.show()
    
# chi2 test to check against a constant profile
def get_chi2(combined_ts, unfolded_ts, folded_ts, period, flux_column_name='flux'):
    # std of unfolded time series - load the data from the different files first, combine and then calculate std
    # x_bar: from Larsson(1996) mean of unfolded data since we want to test against a constant profile (?folded time series)
    # std: from Larsson(1996) std of unfolded data divided by #points in bin = test period / bin width
    
    # mean and std of unfolded data
    mean, std = get_stats_timeseries( combined_ts )
    #mean2, median2, std2 = sigma_clipped_stats(unfolded_ts[flux_column_name])
    print(mean, std)
    n_i = float(get_number_slices(unfolded_ts, period))
    print(f'number of data points in binin chi2 function: {n_i}')
    chi2 = 0
    # total number of bins is given by the length of the folded data
    #print(f'length of folded timeseries:{len(folded_ts)}')
    for i in folded_ts[flux_column_name]:
        # number of data points in bin i - for this case is the same in every bin, but could be different in experimental data if the resolution changes or there are gaps in the data
        #print(f'number of data points in bin: {n_i}')
        #n_i = ((unfolded_ts['time'].value[-1]) - (unfolded_ts['time'].value[0])) / period
        #n_i = float(get_bins_per_slice(unfolded_ts, period))
        #n_i = float(get_number_slices(folded_ts, period))
        #n_i = get_number_slices(folded_ts, period)
        #print(f'number of data points in binin chi2 function: {n_i}')
        chi2 += ( i/n_i - mean )**2 / ( std / np.sqrt(n_i) )**2
        #chi2= 0
    #chi2 = chi2 / (len(folded_ts) - 1)
    return chi2

def get_stats_timeseries(timeseries, column_name='flux'):
    '''calculates mean and std of a timeseries column.'''
    flux = timeseries[column_name]
    return np.mean(flux), np.std(flux)

def plot_timeseries(ts, test_period, filename):
    plt.figure(figsize=(10, 5))
    plt.plot(ts.time.value, ts['flux'], 'o', label='Folded TS', color='blue' )
    plt.xlabel('Time [MJD]')
    plt.ylabel('Amplitude [a.u.]')
    plt.title(f'Folded TimeSeries with Testperiod {test_period}')
    plt.grid(True)
    plt.legend()
    plt.savefig( filename + '.png' )
    plt.show()
    
def load_timeseries(directory=None, filename=None):
    '''returns a list of astropy.timeseries found in directory with the following filename 'timeseries_i.dat', where i is an integer value starting from 0.'''
    if directory is not None:
        file_pattern = os.path.join(directory, 'timeseries_*.dat')
        timeseries_files = glob.glob(file_pattern)
    if filename is not None:
        timeseries_files = glob.glob(filename)
        
    if not timeseries_files:
        print(f"No timeseries files found in the directory: {directory}")
        return

    # Sort files by the integer value i in the filename
    timeseries_files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
    
    if len(timeseries_files) == 1:
        return TimeSeries.read(timeseries_files[0], format='ascii', time_column='time', time_format='mjd')
    else:
        timeseries = []
        for file in timeseries_files:
            ts = TimeSeries.read(file, format='ascii', time_column='time', time_format='mjd')
            timeseries.append(ts)
        return timeseries

def test_testperiod(ts, test_period):
    '''checks if the testperiod is valid and changes it if not.'''
    if Decimal(str(test_period)) >= get_timespan(ts):
        print(f'test period is larger than or equal to the timespan of the timeseries.')
        # period should be max half of the timespan to fold at least two slices but also multiple of binwidth, elseif period > half timespan, fill the timeseries with median
        test_period = ( ( get_timespan(ts) / 2 ) // get_binwidth(ts) ) * get_binwidth(ts)
        print(f'choosing testperiod of {test_period} instead.')
    #if Decimal(str(test_period)) % get_binwidth(ts) == 0:
    elif Decimal(str(test_period)) <= get_binwidth(ts):
        print(f'test period is smaller or equal to the binwidth.')
        #if ( Decimal(str(test_period)) % get_binwidth(ts) ) == Decimal('0.0'): # only true when test_period==binwdith
        test_period = get_binwidth(ts) * 2
        #else:
            #if ( Decimal(str(test_period)) // get_binwidth(ts) + 1 ) == Decimal('1.0'):
            #    test_period = get_binwidth(ts) * 2
            #else:
            #    test_period = ( Decimal(str(test_period)) // get_binwidth(ts) + 1 ) * get_binwidth(ts)
        print(f'choosing testperiod of {test_period} instead.')
    else:
    #elif Decimal(str(test_period)) > ( get_timespan(ts) / 2 ):
        if Decimal(str(test_period)) % get_binwidth(ts) != Decimal('0.0'):
            print(f'test period is not a multiple of the bin width.')
            #test_period = get_binwidth(ts) * 2
        #else:
            #if Decimal(str(test_period)) // get_binwidth(ts) + 1 == Decimal('1.0') :
            # rounding down, to round up: add 1
            test_period = ( Decimal(str(test_period)) // get_binwidth(ts) ) * get_binwidth(ts)
            #else:
            #    test_period = ( Decimal(str(test_period)) // get_binwidth(ts) + 1 ) * get_binwidth(ts)
        #elif get_timespan(ts) % Decimal(str(test_period)) != Decimal('0.0'):
        #elif get_number_bins(ts) % Decimal(str(test_period)) != Decimal('0.0'):
        # the length is the important parameter so the array_split() operation works in the folding function!!
        elif (Decimal(str(len(ts))) % Decimal(str(test_period)) != Decimal('0.0')) or (Decimal(str(len(ts))) % get_bins_per_slice(ts, test_period) != Decimal('0.0')):
            #print(f'test_period is not a submultiple of the timeseries timespan.')
            ### fill the timeseries with mean of entire timeseries
            ts = padd_timeseries(ts, test_period)
            #return Decimal(str(test_period)), padd_timeseries(ts, test_period)
            #print(f'choosing testperiod of {test_period} instead.')
    return Decimal(str(test_period)), ts
            
   # else:
       # if test_period <= get_binwidth(ts):
        #    print(f'test period is smaller or equal to the binwidth. ')
    #if test_period > get_binwidth(ts):
    #    print('test period is valid.')
    #else:
    #    raise ValueError('Test period should be larger than the binwidth of the time series.')
    
def get_stats_timeseries(timeseries, column_name='flux'):
    '''calculates mean and std of a timeseries column.'''
    flux = timeseries[column_name]
    return np.mean(flux), np.std(flux)

def padd_timeseries(ts, test_period, flux_column_name='flux'):
    print(f'padding the timeseries with mean values.')
    #mean, median, std = sigma_clipped_stats(ts)
    mean, std = get_stats_timeseries(ts)
    #add_length = (( get_timespan(ts) // Decimal(str(test_period)) + Decimal('1.0') )* Decimal(str(test_period)) - get_timespan(ts)) / get_binwidth(ts)
    #new_N_bins = get_number_slices(ts, test_period) * Decimal(str(test_period))
    add_length = ((get_number_slices(ts, test_period) * Decimal(str(test_period))) - get_timespan(ts)) / get_binwidth(ts) - 1
    if add_length == 0:
        add_length = 1
    #new_length = ( get_number_slices(ts, test_period) * Decimal(str(test_period)) ) 
    #new_length = diff + get_time
    #print(get_number_slices(ts, test_period), Decimal(str(test_period)), get_timespan(ts), add_length)
    #print( new_length)
    ## remove time column, add new time column, what should N be?
    #add_times = Time(np.linspace(ts.time.value[-1], float(get_binwidth(ts)), int(N), endpoint=False), format='mjd')
    add_ts = TimeSeries(time_start=ts.time[-1] + get_binwidth(ts)*u.d, time_delta=get_binwidth(ts)*u.d, data={flux_column_name: [mean]*int(add_length)})
    #print(add_ts)
    # how to add the flux?? does a fill empty places method exist?
    return vstack([ts, add_ts])

def get_number_bins(ts):
    '''calculates the number of bins in a ts.'''
    return get_timespan(ts) / get_binwidth(ts)

def get_binwidth(ts):
    '''calculates the binwidth of a timeseries.'''
    return Decimal(str(ts.time.value[1])) - Decimal(str(ts.time.value[0]))

def get_timespan(ts):
    '''returns the entire timespan of the timeseries.'''
    return Decimal(str(ts.time.value[-1])) - Decimal(str(ts.time.value[0]))

def get_number_slices(ts, test_period):
    '''returns the integer number of slices to fold depending on the chosen test_period.'''
    
    '''if test_period >= get_timespan(ts):
    #if test_period >= (ts.time.value[-1]) - (reference_time):
        raise ValueError('test period is larger than or equal to the timespan of the timeseries.')
    else:'''
    #print( get_timespan(ts), Decimal(str(test_period)))
    if get_timespan(ts) % Decimal(str(test_period)) != Decimal('0.0'):
         return ( get_timespan(ts) // Decimal(str(test_period)) ) + Decimal('1.0')
    else: 
         return get_timespan(ts) / Decimal(str(test_period))
    

def get_bins_per_slice(ts, test_period):
    #if Decimal(str(test_period)) % get_binwidth(ts) == 0:
    return int( Decimal(str(test_period)) / get_binwidth(ts) )
    #else:
    #    return 

def fold_timeseries(ts, test_period, reference_time=None):
    # set the reference time: the starting time of the folded ts
    if reference_time is None:
        reference_time = ts.time.value[0]
    else:
        reference_time = reference_time
    #print(f'reference_time: {reference_time}')
    #print(len(ts.time), ts.time.value[-1])
    #print( (ts.time.value[1]) - (ts.time.value[0]) )
    binwidth = get_binwidth(ts)
    print(f'binwidth: {binwidth}')
    
    # number of slices
    N = get_number_slices(ts, test_period)
    print( f'Number of slices: {N}' )
    
    test_period, ts = test_testperiod(ts, test_period)
    #print(f'length of ts: {len(ts)}')
    '''
    # both methods work, but introduce rounding errors
    binwidth = get_binwidth(ts)
    print(f'binwidth: {binwidth}')
    
    # number of slices
    N = get_number_slices(ts, test_period)
    print( f'Number of slices: {N}' )
    '''
    
    # number of N=test_period/binwidth bins in one ts slice
    N_bins = get_bins_per_slice(ts, test_period)
    print(f'Number of bins in one slice: {N_bins}')
    
    # set time reference interval
    # introduces rounding errors
    ref_int = Time(np.linspace(0.0, float(test_period), int(N_bins), endpoint=False), format='mjd')
    #print(ref_int, len(ref_int))
    #ref_int = np.linspace(0., test_period, int(N_bins), endpoint=False).tolist()
    #print(ref_int)
    #ref_int = ts.time[]
    #print( ref_int, len(ref_int) )
    
    # split the ts
    #print(ts.as_array())
    splitted_ts = np.array_split(ts.as_array(), int(N))
    #print(splitted_ts)
    #for s in splitted_ts: print(len(s))
    #print(splitted_ts[0], len(splitted_ts[0]))
    #print([x[1] for x in splitted_ts[0]])
    
    # folding the splitted ts
    #print(f'unfolded original timeseries: \n {ts}')
    #'''
    folded_ts = TimeSeries(time=ref_int)
    #folded_ts = TimeSeries(time_start=0.0, time_delta=binwidth*u.d, n_samples=N_bins)
    #print(folded_ts)
    for i, split in enumerate(splitted_ts):
        #print(i, len(ts))
        #print( split[0] )
        #if folded_ts is None:
        if i == 0:
            folded_ts['flux'] = [x[1] for x in split]
            #print(folded_ts)
            #reference_time = folded_ts.time.value[0]
            #folded_ts['time'] = np.round( folded_ts.time.value - reference_time, 2 )
            #print(f'folded timeseries: \n {folded_ts}')
        else:
            folded_ts['flux'] += [x[1] for x in split]
            #folded_ts['flux'] += TimeSeries(data=ts)['flux']
            #folded_ts['time'] -= reference_time * i
        #print(f'folded timeseries: \n {folded_ts}')
    #print(folded_ts[0])
    #plot_timeseries(folded_ts, test_period, 'folded_ts' + str(test_period))
    #'''
    '''
    # fold ts
    folded_ts = TimeSeries(time=times[indices])
    return folded_ts
    '''
    return folded_ts

#'''
def plot_period_slices(ts, periods, filename='bins_per_epoch'):
    plt.figure(figsize=(10, 5))
    y = [get_number_slices(ts, p) for p in periods if p > Decimal('0.0')]
    plt.plot(periods, y, 1, label='Number of bins per epoch', color='blue' )
    plt.xlabel('Time [MJD]')
    plt.ylabel('Amplitude [a.u.]')
    plt.title(f'Number of bins per epoch')
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.savefig( filename + '.png' )
    plt.show()
#'''
if __name__ == '__main__':
    arguments = docopt(__doc__)
    directory = arguments['--directory']
    period = float(arguments['--period'])
    reference_mjd = float(arguments['--reference_mjd'])
   
    #'''
    timeseries = load_timeseries(directory)
    #print(timeseries[0])
    #'''
    #print(f'timespan: {get_timespan(timeseries[0])}')
    '''
    print(f'binwidth: {get_binwidth(timeseries[0])}')
    period, timeseries[0] = test_testperiod(timeseries[0], period)
    print(len(timeseries[0]))
    print(f'test period: {period}')
    print(f'number of slices: {get_number_slices(timeseries[0], period)}')
    print(f'bins per slice: {get_bins_per_slice(timeseries[0], period)}')
    print(f'(padded) timeseries:\n {timeseries[0]}')
    '''
    #folded_ts = fold_timeseries(timeseries[0], period)
    #'''
    
    # combine folded time series into one
    '''use to test one period
    folded_ts = None
    for i, ts in enumerate(timeseries):
        if folded_ts is None:
            folded_ts = fold_timeseries(ts, period, reference_mjd)
            folded_ts['time'] -= reference_mjd
        else:
            folded_ts['flux'] += fold_timeseries(ts, period, reference_mjd)['flux']
        #plot_timeseries(folded_ts, period, 'combined_folded')
    
    plot_timeseries(folded_ts, period, 'combined_folded' + str(period))
    '''
    #ts = combine_timeseries(directory)
    #print(ts)
    
    #'''
    # determine periods from timeseries
    combined_ts = combine_timeseries(directory)
    # the possible periods have the same precision as the binwidth
    # bw < period < timespan of unfolded ts (one file)
    #periods = get_binwidth(combined_ts)*2 + np.arange(Decimal('0.0'), Decimal('200.'))*get_binwidth(combined_ts)
    periods = Decimal('0.0') + np.arange(Decimal('0.0'), Decimal('100.'))*Decimal('0.5')
    print(periods)
    test_periods = periods[1:]
    #plot_period_slices(ts, periods)
    #'''
    chi2s = []
    # calculate chi2 landscape
    with open('chi2_landscape.txt', 'w', newline='') as f:
        #write = csv.writer(f)
        #header = ['period', 'chi2']
        #writer.writerow(header)
        f.write('# period chi2\n')
        for p in test_periods:
            #if get_number_slices()
            
            # for previous test
            #if p >= Decimal('1.99'):
            #    continue
            
            
            # this could be wrong, some test periods are altered and the changes will not appear in the file
            f.write(str(p) + ' ')
            print(f'testperiod: {p}')
            
            folded_ts = None
            for i, ts in enumerate(timeseries):
                #if get_number_slices(ts, p)
                if i == 0:
                    plot_period_slices(ts, test_periods)
                if folded_ts is None:
                    folded_ts = fold_timeseries(ts, p, reference_mjd)
                    folded_ts['time'] -= reference_mjd
                else:
                    folded_ts['flux'] += fold_timeseries(ts, p, reference_mjd)['flux']
            chi2 = get_chi2(combined_ts, ts, folded_ts, p)
            f.write(str(chi2) + '\n')
            f.flush()
            chi2s.append(chi2)
        #plot_timeseries(folded_ts, period, 'combined_folded')
    
            plot_timeseries(folded_ts, p, 'combined_folded' + str(p))
    #chi2 = get_chi2(timeseries[0], folded_ts, period)
    #print(chi2)
    plot_chi2_landscape(chi2s=chi2s, periods=periods, output='chi2_lanscape')
    #'''