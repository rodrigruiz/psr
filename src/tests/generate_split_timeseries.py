#!/usr/bin/env python3
"""
Generate and split time series into multiple files.

Usage:
  script.py --time_start=<time_start> --time_stop=<time_stop> --number_of_points=<number_of_points> --signal_strength=<signal_strength> --period=<period> --num_files=<num_files> --just_noise=<just_noise>
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
plt.style.use('~software/Psr/src/latex.mplstyle')

def gauss(x, mu=0, sigma=1):
    '''returns a normalized gaussian function.'''
    y = 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-(x-mu)**2/2/sigma**2)
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
    #N = int(np.floor((time[-1] - time[0]).to(u.day).value / P)) + 2
    #print(f"Number of pulses N = {N}")
    
    # Index of the phase
    idx_P_0 = int(np.round(P_0 / delta_t.to(u.day).value))
    print(f"idx_P_0={idx_P_0}")
    # Index width of one pulse
    #idx_P = int(np.round(P / delta_t.to(u.day).value))
    #print(f"idx_P={idx_P}")
    # Initialize signal
    PulseTrain = np.sin(2*np.pi/P * t_values + P_0) + 100.
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

# def generate_pulse_train(t, delta_t, P, P_0, sigma, a):
#     # Rescale phase if necessary, so the complete timeseries is filled with pulses
#     # Maps the input P_0 to the interval [-P/2, P/2)
#     P_0 = (P_0 + P/2) % P - P/2

#     # Number of pulses
#     # The floor function yields the number of complete pulses, 
#     # plus the two possible incomplete pulses at the beginning and the end
#     # If no incomplete pulse appears at the beginning or end, the corresponding slice will just yield an empy array
#     # Keep in mind, that a signal with only complete pulses, would still require at leas a '+1' to compensate for floating point errors
#     N = int(np.floor(t[-1]/P)) + 2

#     # Index of the phase
#     idx_P_0 = int(np.round(P_0/delta_t))
    
#     # Indexwidth of one pulse
#     idx_P = int(np.round(P/delta_t))
    
#     # Initialize signal
#     PulseTrain = np.zeros(len(t))

#     # Loop over all Pulses
#     for i in range(0, N):
#         # Calculate lower and upper index for the pulse windows
#         idx_l = idx_P_0 + idx_P//2 + (i-1)*idx_P
#         idx_u = idx_P_0 + idx_P//2 + i*idx_P
#         # Set first index to zero for cut off pulse windows at the beginning
#         if idx_l < 0:
#             idx_l = 0

#         # Calculate pulse peaks for each pulse window
#         PulseTrain[idx_l:idx_u] += a * gauss(t[idx_l:idx_u], i*P + P_0, sigma)

#     return PulseTrain
    

def generate_time_series(time_start, time_stop, number_of_points, signal_strength, period, num_files, just_noise=False):
    times = Time(np.linspace(float(time_start), float(time_stop), int(number_of_points), endpoint=False), format='mjd')

    # delta_t = times[1] - times[0]
    P = float(period)  # Period of the pulse train
    P_0 = 0 # Phase offset
    
    a = float(signal_strength)  # Amplitude
    
    # gaussian signal
    sigma = 1 #Standard deviation
    flux = generate_pulse_train_gauss(times, P, P_0, sigma, a) 
    label = f'Gaussian Pulse Train' 
    title = f'Gaussian Pulse Train with Period of {P} and Sigma of {sigma}'
    
    # sinusodial signal
    #flux = generate_pulse_train_sin(times, P, P_0, a)
    #label = 'Sinusodial Pulse Train'
    #title = f'Sinusodial Pulse Train with Period of {P}'
    
    # for visualization - can be removed later
    plt.figure(figsize=(10, 5))
    plt.plot(times.value, flux, label=label, color='blue' )
    plt.xlabel('Time [MJD]')
    plt.ylabel('Amplitude [a.u.]')
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.savefig( 'pulse_train.png' )
    plt.show()
    
    # std defines noise strength
    # poisson
    noise = np.random.poisson(100, int(number_of_points))
    title_poisson = ' and Poisson distributed Noise'
    
    # gaussian
    #noise = np.random.normal(0, 0.5, int(number_of_points))
    #title_gauss = ' and gaussian distributed Noise'
    
    if just_noise:
        flux = noise
    else:
        flux += noise
    
    # for visualization - can be removed later
    plt.figure(figsize=(10, 5))
    plt.plot(times.value, flux, label=label + ' + Noise', color='blue' )
    plt.xlabel('Time [MJD]')
    plt.ylabel('Amplitude [a.u.]')
    plt.title(title + title_poisson)
    plt.grid(True)
    plt.legend()
    plt.savefig( 'pulse_train+noise.png' )
    plt.show()
    
    # for visualization - can be removed later
    plt.figure(figsize=(10, 5))
    plt.plot(times.value[0:250], flux[0:250], label=label + ' + Noise', color='blue' )
    plt.xlabel('Time [MJD]')
    plt.ylabel('Amplitude [a.u.]')
    plt.title(title + title_poisson)
    plt.grid(True)
    plt.legend()
    plt.savefig( 'pulse_train+noise_zoom.png' )
    plt.show()

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
    plt.savefig(f'timeseries.png')
    plt.show()
    

if __name__ == '__main__':
    arguments = docopt(__doc__)
    generate_time_series(
        arguments['--time_start'],
        arguments['--time_stop'],
        arguments['--number_of_points'],
        arguments['--signal_strength'],
        arguments['--period'],
        arguments['--num_files'],
        arguments['--just_noise']
    )
