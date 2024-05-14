import numpy as np
import h5py
from astropy.timeseries import TimeSeries
from astropy.time import Time
import astropy.units as u
from docopt import docopt
from plens.TimeSeries import saveTimeSeries
import matplotlib.pyplot as plt
#plt.style.use('~/software/Psr/src/latex.mplstyle')

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

def plot_pulsetrain(times, counts, output='pulsetrain', color='blue', label='PulseTrain', title='PulseTrain'):
    plt.figure(figsize=(10, 5))
    plt.plot(times, counts, label=label, color=color )
    plt.xlabel('Time [MJD]')
    plt.ylabel('Counts [a.u.]')
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.savefig( output + '.png' )
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
        title = r'Gaussian Pulse Train with Frequency $f={}$ and Sigma $\sigma={}$'.format(f, sigma)
        output = 'pulsetrain_gauss'
    elif pulseshape == 'mvm':
        kappa = 10.
        flux = generate_pulse_train_peak(times, P, P_0, pulseshape, a/8, kappa=kappa)
        label = f'modified von Mises Pulse Train' 
        title = r'modified von Mises Pulse Train with Frequency $f={}$ and $\kappa={}$'.format(f, kappa)
        output = 'pulsetrain_mvm'
    elif pulseshape == 'sin':
        flux = generate_pulse_train_sin(times, P, P_0, a, 50.)
        label = f'Sinusodial Pulse Train' 
        title = r'Sinusodial Pulse Train with Frequency $f={}$'.format(P, sigma)
        output = 'pulsetrain_sin'
    else:
        raise ValueError(
            "pulseshape should be `gauss`, `mvm` or `sin`!"
        )
        
    # plot the pure signal
    plot_pulsetrain(times.value, flux, output=output, label=label, title=title)
    
    # std defines noise strength
    # poisson noise
    if noise == 'poisson':
        noise = np.random.poisson(50, int(number_of_points))
        label_noise = ' + poissonian distributed Noise'
        output_noise = '+poissonnoise'
    elif noise == 'gauss':
        noise = np.random.normal(0, 0.5, int(number_of_points))
        label_noise = ' + gaussian distributed Noise'
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
        flux += noise
    #flux = noise
    
    # plot the pulsetrain with noise
    plot_pulsetrain(times.value, flux, output=output+output_noise, label=label+label_noise, title=title+label_noise)
    # zoom into the pulsetrain
    plot_pulsetrain(times.value[:250], flux[:250], output=output+output_noise+'_zoom', label=label+label_noise, title=title+label_noise)
    
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
        
