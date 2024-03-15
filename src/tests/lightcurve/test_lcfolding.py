# %load_ext autoreload
# %autoreload 2
# %matplotlib notebook

import numpy as np
from stingray.pulse.search import epoch_folding_search, z_n_search
import matplotlib.pyplot as plt
import seaborn as sb
import matplotlib as mpl
from stingray.pulse.pulsar import fold_events
from stingray.pulse.search import plot_profile
from stingray import Lightcurve
mpl.rcParams['figure.figsize'] = (10, 6)
from astropy.time import Time
from stingray import EventList, Lightcurve
import astropy.units as u

def sinusoid(times, frequency, baseline, amplitude, phase):
    return baseline + amplitude * np.sin(2 * np.pi * (frequency * times + phase))

def generate_pulse_train_peak(time, f, P_0, pulseshape, a, baseline, sigma=None, kappa=None):
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
    
    #t_values = time
    
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
            PulseTrain[idx_l:idx_u] += MVMD(t_values[idx_l:idx_u], f, baseline, a, P_0, kappa)
        else:
            raise ValueError(
                "pulseshape can only be `gauss` or `mvm`!"
            )
            
    return PulseTrain

def MVMD(t, f, baseline, a, phi, kappa):
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
    
    y = a * ( np.exp(kappa*np.cos(2*np.pi*f*t+phi))-np.exp(-kappa))/(np.i0(kappa) - np.exp(-kappa)) + baseline
    return y

period = 1.203501
mean_countrate = 50
pulsed_fraction = 1.
bin_time = 0.1
obs_length = 3000

t0 = np.random.uniform( 1489328336, 1489328336+obs_length, 1000 )
print(len(t0))
t = np.arange(1489328336, 1489328336+obs_length, bin_time)
#print(t[0])

# The continuous light curve
#counts = sinusoid(t, 1 / period, mean_countrate, 
#                  0.5 * mean_countrate * pulsed_fraction, 0) * bin_time
counts = MVMD(t, 1 / period, mean_countrate, 0.5 * mean_countrate * pulsed_fraction, np.pi, 5.)
len(counts)
#counts = generate_pulse_train_peak(Time(t, format='unix'), 1 / period, 0, 'mvm', 0.5 * mean_countrate, mean_countrate, kappa=5.)
lc = Lightcurve(t, counts, dt=bin_time, skip_checks=True)

ev = EventList()
ev.simulate_times(lc)
print(len(ev.time), np.shape(ev.time)) #np.array

new = np.sort(np.concatenate((t0, ev.time)))
print(len(new), np.shape(new))

plt.figure()
plt.plot(t[:500], counts[:500], color='blue' )
plt.xlabel('Time [a.u.]')
plt.ylabel('Counts [a.u.]')
plt.title('Sinusoid')
plt.grid(True)
#plt.legend()
plt.savefig( 'MVM_stingray.png' )
plt.show()

'''
nbin = 32

ph, profile, profile_err = fold_events(events.time, 1/period, nbin=nbin)
_ = plot_profile(ph, profile)

ph, profile, profile_err = fold_events(events.time, 1/1.1, nbin=nbin)
_ = plot_profile(ph, profile)
'''