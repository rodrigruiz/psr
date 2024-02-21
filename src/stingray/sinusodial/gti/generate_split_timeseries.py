#!/usr/bin/env python3
"""
Generate and split time series into multiple files.

Usage:
  script.py --time_start=<time_start> --time_stop=<time_stop> --number_of_points=<number_of_points> --signal_strength=<signal_strength> --period=<period> --num_files=<num_files>
  script.py (-h | --help)

Options:
  -h --help               Show this help message and exit.
  --time_start=<time_start>     Start time
  --time_stop=<time_stop>       Stop time
  --number_of_points=<number_of_points> Number of points
  --signal_strength=<signal_strength>   Signal strength
  --period=<period>           Period of the sinusoid signal
  --num_files=<num_files>     Number of output files
"""

import numpy as np
from astropy.timeseries import TimeSeries
from astropy.time import Time
import astropy.units as u
from docopt import docopt
import matplotlib.pyplot as plt
#plt.style.use('../latex.mplstyle')

def gauss(x, mu=0, sigma=1):
    '''returns a normalized gaussian function.'''
    y = 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-(x-mu)**2/2/sigma**2)
    return y

def sinusoid(times, period, baseline, amplitude, phase):
    return baseline + amplitude * np.sin(2 * np.pi * (times / period + phase))

def MVMD(t, f, phi, kappa, a):
    """Modivied von Mises distribution (MVMD).
        Reference: "Fourier Techniques for Very Long Astrophysical Time-Series Analysis" 
                    Scott M. Ransom et al 2002 AJ 124 1788
                    
    Parameters
    ----------
        t : array_like
            Values to evaluate the MVMD at.
        kappa : float
            Shape parameter giving the width of the function.
        a : float
            Amplitude of the Pulse. Equates to the area of one pulse.
        f : float
            Frequency of the pulse train.
        phi : float
            Phase offset of the pulse train.
    
    Returns
    -------
        y : 1darray
            Evaluated values of the MVMD.
    
    For kappa -> 0 the MVMD converges towards a sinusod.
    For kappa -> infinity the MVMD converges towards a Gaussian. 1/kappa corresponds to sigma^2. Keep in mind, that  
        
    """
    
    y = a * ( np.exp(kappa*np.cos(2*np.pi*f*t+phi))-np.exp(-kappa))/(np.i0(kappa) - np.exp(-kappa))
    return y

def generate_pulse_train_gauss(time, P, P_0, sigma, a):
    '''generates a pulse train in the entire 'time' interval with period P.
each pulse is represented by a gaussian function with amplitude a and sigma.
returns: np.array containing pulse train'''
    delta_t = time[1] - time[0]
    # Extract values from the Time object
    t_values = time.value
    # Rescale phase if necessary, so the complete timeseries is filled with pulses: # Maps the input P_0 to the interval [-P/2, P/2)
    P_0 = (P_0 + P/2) % P - P/2
    print(f'Phase Offset P_0 = {P_0}')

    # Number of pulses
    # The floor of the scalar x is the largest integer i, such that i <= x
    N = int(np.floor((time[-1] - time[0]).to(u.day).value / P)) + 2
    print(f"Number of pulses N = {N}")
    
    # Index of the phase
    idx_P_0 = int(np.round(P_0 / delta_t.to(u.day).value))
    print(f"idx_P_0={idx_P_0}")
    # Index width of one pulse
    idx_P = int(np.round(P / delta_t.to(u.day).value))
    print(f"idx_P={idx_P}")
    # Initialize signal
    PulseTrain = np.zeros(len(time))
    # Loop over all Pulses
    for i in range(0, N):
        # Calculate lower and upper index for the pulse windows
        idx_l = idx_P_0 + idx_P // 2 + (i - 1) * idx_P
        idx_u = idx_P_0 + idx_P // 2 + i * idx_P
        #print(f"idx_l={idx_l}  idx_u={idx_u}")
        # Set first index to zero for cut off pulse windows at the beginning
        if idx_l < 0:
            idx_l = 0
        # Calculate pulse peaks for each pulse window
        PulseTrain[idx_l:idx_u] += a * gauss(t_values[idx_l:idx_u], i * P + P_0 + t_values[0], sigma)
    return PulseTrain

def generate_pulse_train_sin(time, P, P_0, a):
    '''generates a pulse train in the entire 'time' interval with period P.
each pulse is represented by a sine function with amplitude a and sigma.
returns: np.array containing pulse train'''
    delta_t = time[1] - time[0]
    # Extract values from the Time object
    t_values = time.value
    # Rescale phase if necessary, so the complete timeseries is filled with pulses: # Maps the input P_0 to the interval [-P/2, P/2)
    P_0 = (P_0 + P/2) % P - P/2
    print(f'Phase Offset P_0 = {P_0}')

    # Number of pulses
    # The floor of the scalar x is the largest integer i, such that i <= x
    #N = int(np.floor((time[-1] - time[0]).to(u.day).value / P)) + 2
    #print(f"Number of pulses N = {N}")
    
    # Index of the phase
    idx_P_0 = int(np.round(P_0 / delta_t.to(u.day).value))
    print(f"idx_P_0={idx_P_0}")
    # Index width of one pulse
    #idx_P = int(np.round(P / delta_t.to(u.day).value))
    #print(f"idx_P={idx_P}")
    # Initialize signal
    PulseTrain = a * np.sin(2*np.pi/P * t_values + P_0) + 50.
    #multiple sine
    PulseTrain += a/2 * np.sin(2*np.pi / (P/3) * t_values + P_0) + 50.
    
    # for visualization - can be removed later
    plt.figure(figsize=(10, 5))
    plt.plot(t_values[:200], PulseTrain[:200], color='blue' )
    plt.plot(t_values[:200], a * np.sin(2*np.pi/P * t_values[:200] + P_0) + 50., color='red')
    plt.plot(t_values[:200], a/2 * np.sin(2*np.pi / (P/3) * t_values[:200] + P_0) + 50., color='red', ls='--')
    plt.xlabel('Time [MJD]')
    plt.ylabel('Amplitude [a.u.]')
    #plt.title(title)
    plt.grid(True)
    #plt.legend()
    plt.savefig( 'pure_pulse_2sine.png' )
    plt.show()
    '''
    # using stingray function
    #PulseTrain = sinusoid(t_values, P, 100., a, P_0)
    mean_countrate = 50.
    pulsed_fraction = 0.2 
    bin_time = 0.01
    PulseTrain = sinusoid(t_values, P, mean_countrate, 
                  0.5 * mean_countrate * pulsed_fraction, P_0) * bin_time
    '''
    '''
    # Loop over all Pulses
    for i in range(0, N):
        # Calculate lower and upper index for the pulse windows
        idx_l = idx_P_0 + idx_P // 2 + (i - 1) * idx_P
        idx_u = idx_P_0 + idx_P // 2 + i * idx_P
        #print(f"idx_l={idx_l}  idx_u={idx_u}")
        # Set first index to zero for cut off pulse windows at the beginning
        if idx_l < 0:
            idx_l = 0
        # Calculate pulse peaks for each pulse window
        PulseTrain[idx_l:idx_u] += a * gauss(t_values[idx_l:idx_u], i * P + P_0 + t_values[0], sigma)
    '''
    return PulseTrain


def generate_pulse_train_mvm(time, P, P_0, a):
    '''generates a pulse train in the entire 'time' interval with period P.
each pulse is represented by a gaussian function with amplitude a and sigma.
returns: np.array containing pulse train'''
    delta_t = time[1] - time[0]
    # Extract values from the Time object
    t_values = time.value
    # Rescale phase if necessary, so the complete timeseries is filled with pulses: # Maps the input P_0 to the interval [-P/2, P/2)
    P_0 = (P_0 + P/2) % P - P/2
    print(f'Phase Offset P_0 = {P_0}')

    # Number of pulses
    # The floor of the scalar x is the largest integer i, such that i <= x
    N = int(np.floor((time[-1] - time[0]).to(u.day).value / P)) + 2
    print(f"Number of pulses N = {N}")
    
    # Index of the phase
    idx_P_0 = int(np.round(P_0 / delta_t.to(u.day).value))
    print(f"idx_P_0={idx_P_0}")
    # Index width of one pulse
    idx_P = int(np.round(P / delta_t.to(u.day).value))
    print(f"idx_P={idx_P}")
    # Initialize signal
    PulseTrain = np.zeros(len(time))
    # Loop over all Pulses
    for i in range(0, N):
        # Calculate lower and upper index for the pulse windows
        idx_l = idx_P_0 + idx_P // 2 + (i - 1) * idx_P
        idx_u = idx_P_0 + idx_P // 2 + i * idx_P
        #print(f"idx_l={idx_l}  idx_u={idx_u}")
        # Set first index to zero for cut off pulse windows at the beginning
        if idx_l < 0:
            idx_l = 0
        # Calculate pulse peaks for each pulse window
        #PulseTrain[idx_l:idx_u] += a * gauss(t_values[idx_l:idx_u], i * P + P_0 + t_values[0], sigma)
        PulseTrain[idx_l:idx_u] += MVMD(t_values[idx_l:idx_u], 1/P, P_0, 10, a)
    '''
    if gti is None:
            PulseTrain[idx_l:idx_u] += MVMD(t_values[idx_l:idx_u], 1/P, P_0, 10, a)
        else:
            for t in t_values[idx_l:idx_u]:
                if gti[0] <= t <= gti[1]:
                    continue
                else:
                    PulseTrain[idx_l:idx_u] += MVMD(t_values[idx_l:idx_u], 1/P, P_0, 10, a)
    '''
    #for t, i in enumerate(t_values):
        
        # for visualization - can be removed later
    plt.figure(figsize=(10, 5))
    plt.plot(t_values[:200], PulseTrain[:200], color='blue' )
    #plt.plot(t_values[:200], a * np.sin(2*np.pi/P * t_values[:200] + P_0) + 50., color='red')
    #plt.plot(t_values[:200], a/2 * np.sin(2*np.pi / (P/3) * t_values[:200] + P_0) + 50., color='red', ls='--')
    plt.xlabel('Time [MJD]')
    plt.ylabel('Amplitude [a.u.]')
    #plt.title(title)
    plt.grid(True)
    #plt.legend()
    plt.savefig( 'pure_pulse_mvm.png' )
    plt.show()
    
    return PulseTrain
    

def generate_time_series(time_start, time_stop, number_of_points, signal_strength, period, num_files):
    times = Time(np.linspace(float(time_start), float(time_stop), int(number_of_points), endpoint=False), format='mjd')

    # delta_t = times[1] - times[0]
    P = float(period)  # Period of the pulse train
    P_0 = 0 # Phase offset
    
    a = float(signal_strength)  # Amplitude
    
    # gaussian signal
    #sigma = 1 #Standard deviation
    #flux = generate_pulse_train_gauss(times, P, P_0, sigma, a) 
    #label = f'Gaussian Pulse Train' 
    #title = f'Gaussian Pulse Train with Period of {P} and Sigma of {sigma}'
    
    # sinusodial signal
    #flux = generate_pulse_train_sin(times, P, P_0, a)
    #label = 'Sinusodial Pulse Train'
    #title = f'Sinusodial Pulse Train with Period of {P}'
    
    #MVM
    #flux = generate_pulse_train_mvm(times, P, P_0, a)
    flux = 50.
    label = 'modified von Mises Pulse Train'
    title = f'modified von Mises Pulse Train with period {P}'
    
    '''
    # for visualization - can be removed later
    plt.figure(figsize=(10, 5))
    plt.plot(times.value, flux, label=label, color='blue' )
    plt.xlabel('Time [MJD]')
    plt.ylabel('Amplitude [a.u.]')
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.savefig( 'pulse_train_mvm.png' )
    plt.show()
    '''
    # std defines noise strength
    # poisson
    noise = np.random.poisson(50, int(number_of_points))
    title_poisson = ' and Poisson distributed Noise'
    
    # gaussian
    #noise = np.random.normal(0, 0.5, int(number_of_points))
    #title_gauss = ' and gaussian distributed Noise'
    
    
    flux += noise
    '''
    # for visualization - can be removed later
    plt.figure(figsize=(10, 5))
    plt.plot(times.value, flux, label=label + ' + Noise', color='blue' )
    plt.xlabel('Time [MJD]')
    plt.ylabel('Amplitude [a.u.]')
    plt.title(title + title_poisson)
    plt.grid(True)
    plt.legend()
    plt.savefig( 'pulse_train+noise_mvm.png' )
    plt.show()
    
    # for visualization - can be removed later
    plt.figure(figsize=(10, 5))
    plt.plot(times.value[0:250], flux[0:250], label=label + ' + Noise', color='blue' )
    plt.xlabel('Time [MJD]')
    plt.ylabel('Amplitude [a.u.]')
    plt.title(title + title_poisson)
    plt.grid(True)
    plt.legend()
    plt.savefig( 'pulse_train+noise_zoom_mvm.png' )
    plt.show()
    '''
    # Split the time series into multiple files
    file_indices = np.array_split(np.arange(int(number_of_points)), int(num_files))
    plt.figure(figsize=(10, 6))
    
    for i, indices in enumerate(file_indices):
        ts = TimeSeries(time=times[indices])
        ts['counts'] = flux[indices]
        print(f'timeseries: {ts}')
        plt.plot(ts.time.mjd, ts['counts'], '.', markersize=1, label=f'Time Series {i}')
        file_name = f'timeseries_{i}.dat'
        ts.write(file_name, format='ascii.ecsv', overwrite=True)
    plt.xlabel('Time [MJD]')
    plt.ylabel('Counts [a.u.]')
    plt.title(f'period = {period} days.')
    plt.savefig(f'timeseries_mvm.png')
    plt.show()
    

if __name__ == '__main__':
    arguments = docopt(__doc__)
    generate_time_series(
        arguments['--time_start'],
        arguments['--time_stop'],
        arguments['--number_of_points'],
        arguments['--signal_strength'],
        arguments['--period'],
        arguments['--num_files']
    )
