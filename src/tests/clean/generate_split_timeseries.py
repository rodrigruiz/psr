#!/usr/bin/env python3
"""
Generate and split time series into multiple files.

Usage:
  script.py --time_start=<time_start> --time_stop=<time_stop> --number_of_points=<number_of_points> --signal_strength=<signal_strength> --frequency=<frequency> --num_files=<num_files> --pulseshape=<pulseshape> --noise=<noise> [--just_noise]
  script.py (-h | --help)

Options:
  -h --help               Show this help message and exit.
  --time_start=<time_start>     Start time
  --time_stop=<time_stop>       Stop time
  --number_of_points=<number_of_points> Number of points
  --signal_strength=<signal_strength>   Signal strength
  --frequency=<frequency>           Frequency of the sinusoid signal
  --num_files=<num_files>     Number of output files
  --pulseshape=<pulseshape>   Shape of pulses in the pulse train
  --noise=<noise>             Distribution of noise added to pulse train
  --just_noise   Select whether the signal is just noise. [default: False]
"""

import numpy as np
import h5py
from astropy.timeseries import TimeSeries
from astropy.time import Time
import astropy.units as u
from docopt import docopt
from plens.TimeSeries import saveTimeSeries
import matplotlib.pyplot as plt
plt.style.use('~/software/Psr/src/latex.mplstyle')
from plens.plot_latex_size import set_size
width = 418.25368

def gauss(x, mu=0, sigma=1):
    """normalized gaussian function.
    
    Parameters
    ----------
        x : array_like or scalar
            Values to evaluate gauss function at.
        mu : float
            Mean value of gauss function.
            Defaults to 0.
            
        sigma : float
            Standard deviation of gauss function.
            Defaults to 1.
            
    Returns
    -------
        1darray or scalar
        Gauss function.
        
    """
    
    return 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-(x-mu)**2/2/sigma**2)

def sinusoid(times, frequency, baseline, amplitude, phase):
    return baseline + amplitude * np.sin(2 * np.pi * (times * frequency + phase))

def MVMD(t, f, phi, kappa, a):
    """Modified von Mises distribution (MVMD).
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

def generate_pulse_train_peak(time, f, P_0, pulseshape, a, sigma=None, kappa=None):
    """generates a pulse train in the entire `time` interval with pulse frequency `f`.

    Parameters
    ----------
        time : array_like, astropy.time.Time
            Time bins of the pulse train.
            
        f : float
            Frequency of the pulse train.
            
        P_0 : float
            Phase of the pulses.
        
        pulseshape : string
            If 'gauss', each pulse is represented by a gaussian function with amplitude `a` and standard deviation `sigma`.
            If 'mvm', each pulse is represented by a modified von Mises distribution with amplitude `a` and shape parameter `kappa`.
        
        a : float
            Amplitude of pulses.
        
        sigma : float
            Standard deviation of gauss function.
            To be set, if `pulseshape` is 'gauss'.
            
        kappa : float
            Shape parameter of modified von Mises distribution.
            To be set if `pulseshape` is 'mvm'.
        
    Returns
    -------
        np.array 
        Pulse train
        
    """
    
    # binwidth 
    delta_t = time[1] - time[0]
    
    # Extract values from the Time object
    t_values = time.value
    
    # Rescale phase if necessary, so the complete timeseries is filled with pulses: # Maps the input P_0 to the interval [-P/2, P/2)
    #P_0 = (P_0 + P/2) % P - P/2
    #print(f'Phase Offset P_0 = {P_0}')

    # Number of pulses
    # The floor of the scalar x is the largest integer i, such that i <= x
    N = int(np.floor((time[-1] - time[0]).to(u.day).value * f)) + 2
    print(f"Number of pulses N = {N}")
    
    # Index of the phase
    idx_P_0 = int(np.round(P_0 / delta_t.to(u.day).value))
    print(f"idx_P_0={idx_P_0}")
    
    # Index width of one pulse
    idx_P = int(np.round(1 / f / delta_t.to(u.day).value))
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
        if pulseshape == 'gauss':
            PulseTrain[idx_l:idx_u] += a * gauss(t_values[idx_l:idx_u], i * 1/f + P_0 + t_values[0], sigma)
        elif pulseshape == 'mvm':
            PulseTrain[idx_l:idx_u] += MVMD(t_values[idx_l:idx_u], f, P_0, kappa, a)
        else:
            raise ValueError(
                "pulseshape can only be `gauss` or `mvm`!"
            )
            
    return PulseTrain

def generate_pulse_train_sin(time, f, P_0, a, baseline):
    """generates a pulse train in the entire `time` interval with pulse frequency `f`.
   
    Parameters
    ----------
        time : array_like, astropy.time.Time
            Time bins of the pulse train.
            
        f : float
            Frequency of the pulse train. In this case equals the frequency of the sine function.
            
        P_0 : float
            Phase of the pulses.
        
        a : float
            Amplitude of pulses.
            
        baseline : float
            Offset of the pulse train. Shifts the signal on the y-axis.
        
        
    Returns
    -------
        np.array 
        Pulse train
        
    """
    
    # bindwidth
    delta_t = time[1] - time[0]
    
    # Extract values from the Time object
    t_values = time.value
    
    # Rescale phase if necessary, so the complete timeseries is filled with pulses: # Maps the input P_0 to the interval [-P/2, P/2)
    P_0 = (P_0 + P/2) % P - P/2
    print(f'Phase Offset P_0 = {P_0}')
    
    # Initialize signal
    PulseTrain = a * np.sin(2*np.pi*f * t_values + P_0) + baseline
    
    #multiple sine
    #PulseTrain += a/2 * np.sin(2*np.pi / (P/3) * t_values + P_0) + 50.
    
    return PulseTrain

def plot_pulsetrain(times, counts, output='pulsetrain', color='blue', label=None, title='PulseTrain', fraction=1):
    #plt.figure(figsize=(5,3))
    plt.figure(figsize=set_size(width,fraction=fraction))
    plt.plot(times, counts, label=label, color=color )
    plt.xlabel('Time $t$ [s]')
    #plt.xlabel('Time [MJD]')
    plt.ylabel('Counts [a.u.]')
    plt.title(title)
    plt.grid(True)
    if label is not None:
        plt.legend()
    plt.savefig( output + '.pdf' )
    plt.show()

def generate_time_series(time_start, time_stop, number_of_points, signal_strength, frequency, num_files, pulseshape, noise, just_noise=False):
    """TimeSeries (or binned lightcurve).
    
    Parameters
    ----------
        time_start : float
            Starttime of the TimeSeries in units of modified julian date (mjd).
            
        time_stop : float
            Stoptime of the TimeSeries in units of modified julian date (mjd).
            
        number_of_points : int
            Number of points between `time_start` and `time_stop`. 
            Defines the binwidth of the TimeSeries (= time_start - time_stop / number_of_points)
            
        signal_strength : float
            Amplitude of the pulses in the pulse train.
            
        frequency : float
            Frequencyo of the pulses in the pulse train.
            
        num_files : int
            Number of files to split the entire TimeSeries into.
            
        pulseshape : string
            If 'gauss', the signal will be a gaussina pulse train.
            If 'mvm', the signal will be a pulse train with modified von Mises pulses.
            If 'sin', the signal will be a sine function.
            
        noise : string
            If 'poisson', poisson distributed noise will be added to the pulse train.
            If 'gauss', white noise will be added to the pulse train.
        
    Returns
    -------
        Saves the (splitted) astropy.timeseries.TimeSeries to files (timeseries_*.dat) in ascii.ecsv format.
        
    """
    
    times = Time(np.linspace(float(time_start), float(time_stop), int(number_of_points), endpoint=False), format='mjd')

    f = float(frequency)  # Frequency of the pulse train
    P_0 = 0 # Phase offset
    
    a = float(signal_strength)  # Amplitude
    
    if pulseshape == 'gauss':
        sigma = 1. #Standard deviation
        flux = generate_pulse_train_peak(times, f, P_0, pulseshape, a, sigma=sigma)
        label = f'Gaussian Pulse Train' 
        title = r'Gaussian Pulse Train with frequency $f={} \,\mathrm{Hz}$ and Sigma $\sigma={}$'.format(f, sigma)
        output = 'pulsetrain_gauss'
    elif pulseshape == 'mvm':
        kappa = 5.
        #flux = generate_pulse_train_peak(times, f, P_0, pulseshape, a, kappa=kappa)
        flux = MVMD(times.value, f, P_0, kappa, a)
        label = f'modified von Mises Pulse Train' 
        #label = ''
        #title = r'modified von Mises pulsetrain \newline with frequency $f={} \,\mathrm{{Hz}}$ and $\kappa={}$'.format(f, kappa)
        title = r'$f_\mathrm{{MVMD}}(t; f = {} \,\mathrm{{Hz}}, \kappa = {})$'.format(f, kappa)
        title_smaller = r'$f_\mathrm{{MVMD}}(t; f, \kappa)$'
        output = 'pulsetrain_mvm'
    elif pulseshape == 'sin':
        flux = generate_pulse_train_sin(times, f, P_0, a, 50.)
        label = f'Sinusodial Pulse Train' 
        title = r'Sinusodial Pulse Train with Frequency $f={}$'.format(P, sigma)
        output = 'pulsetrain_sin'
    else:
        raise ValueError(
            "pulseshape should be `gauss`, `mvm` or `sin`!"
        )
    
    #print(flux)
    # plot the pure signal, for using hertz -58000. and s on the x-axis
    plot_pulsetrain(times.value - 58000., flux, output=output, title=title, fraction=0.5)
    #plot_pulsetrain(times.value - 58000., flux, output=output, label=label, title=title)
    # std defines noise strength
    # poisson noise
    if noise == 'poisson':
        noise = np.random.poisson(max(flux)*3, int(number_of_points))
        #label_noise = '\n and Poisson noise'
        label_noise = ' and Poisson noise'
        output_noise = '+poissonnoise'
    elif noise == 'gauss':
        noise = np.random.normal(0, 0.5, int(number_of_points))
        label_noise = ' + Gaussian distributed Noise'
        output_noise = '+gaussnoise'
    elif noise == 'exp':
        l = np.log(2) / len(times.value) / 2
        noise = ( 7 * 10 * np.exp(- (times.value - times.value[0]) * l) + np.random.poisson(20, int(number_of_points)) ) #/ 10000 #- (7 * 10**3 * 0.9)
        label_noise = ' + exponentially decaying poisson Noise'
        output_noise = '+exp+poissonnoise'
    else:
        raise ValueError(
            "noise should be either `poisson` or `gauss`!"
        )
    
    # add noise to pulse train
    if just_noise:
       # print(type(noise))
        #flux = [x+1 for x in noise]
        if gauss:
            flux = np.add(noise, 10)
        else:
            flux = noise
    else:
        #print('I should be here')
        flux += noise
    #flux = noise
    
    label=label+label_noise
    # plot the pulsetrain with noise
    plot_pulsetrain(times.value -58000., flux, output=output+output_noise, title=title_smaller+label_noise, fraction=0.5)
    # zoom into the pulsetrain
    plot_pulsetrain(times.value[:250]-58000., flux[:250], output=output+output_noise+'_zoom', title=title+label_noise)
    
    # Split the time series into multiple files
    file_indices = np.array_split(np.arange(int(number_of_points)), int(num_files))
    
    for i, indices in enumerate(file_indices):
        ts = TimeSeries(time=times[indices])
        ts['counts'] = flux[indices]
        
        file_name = f'timeseries_{i}.dat'
        ts.write(file_name, format='ascii.ecsv', overwrite=True)

def split_timeseries(ts, number_of_points, num_files, filename, format='ascii.ecsv', antares=False, timeslice_duration=None):
    file_indices = np.array_split(np.arange(int(number_of_points)), int(num_files))
            
    for i, indices in enumerate(file_indices):
        
        if format == 'hdf5':
            output = filename + f'_{i}.hdf5'
        else:
            output = filename + f'_{i}.dat'
        
        if antares:
            with h5py.File(output, 'w') as output_file:
                saveTimeSeries(ts[indices], timeslice_duration, output_file)
        else:
            ts[indices].write(output, format=format, overwrite=True)
        
if __name__ == '__main__':
    arguments = docopt(__doc__)
    print(arguments)
    generate_time_series(
        arguments['--time_start'],
        arguments['--time_stop'],
        arguments['--number_of_points'],
        arguments['--signal_strength'],
        arguments['--frequency'],
        arguments['--num_files'],
        arguments['--pulseshape'],
        arguments['--noise'],
        arguments['--just_noise']
    )
