# %load_ext autoreload
# %autoreload 2
# %matplotlib notebook

import numpy as np
from stingray.pulse.search import epoch_folding_search, z_n_search
import matplotlib.pyplot as plt
import seaborn as sb
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (10, 6)

#from generate_split_timeseries import generate_pulse_train_gauss, gauss
from stingray import Lightcurve
from stingray.events import EventList
from stingray.pulse.pulsar import fold_events
from stingray.stats import fold_detection_level
from stingray.pulse.search import plot_profile, search_best_peaks
from stingray.pulse.modeling import fit_gaussian, fit_sinc

def sinusoid(times, frequency, baseline, amplitude, phase):
    return baseline + amplitude * np.sin(2 * np.pi * (frequency * times + phase))

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
    t_values = time
    # Rescale phase if necessary, so the complete timeseries is filled with pulses: # Maps the input P_0 to the interval [-P/2, P/2)
    P_0 = (P_0 + P/2) % P - P/2
    print(f'Phase Offset P_0 = {P_0}')

    # Number of pulses
    # The floor of the scalar x is the largest integer i, such that i <= x
    N = int(np.floor((time[-1] - time[0]) / P)) + 2
    print(f"Number of pulses N = {N}")
    
    # Index of the phase
    idx_P_0 = int(np.round(P_0 / delta_t))
    print(f"idx_P_0={idx_P_0}")
    # Index width of one pulse
    idx_P = int(np.round(P / delta_t))
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

'''
period = 1.203501
mean_countrate = 50
pulsed_fraction = 0.2
bin_time = 0.01
obs_length = 3000
''' 
period = 10. # s
P_0 = 0
mean_countrate = 50
#pulsed_fraction = 0.2
sigma = 0.03
bin_time = 0.1
obs_length = 1000

'''
t = np.arange(0, obs_length, bin_time)

# The continuous light curve
#counts = sinusoid(t, 1 / period, mean_countrate, 0.5 * mean_countrate * pulsed_fraction, 0) * bin_time
counts = generate_pulse_train_gauss(t, period, P_0, sigma, mean_countrate )
print(counts)

lc = Lightcurve(t, counts, gti=[[-bin_time / 2, obs_length + bin_time / 2]],
                dt=bin_time, mjdref=t[0], skip_checks=True)

# use the light curve above to simulate an event list for this pulsar.
events = EventList()
events.simulate_times(lc)
print(events.time)
'''

events = EventList().read('./gti/combined_eventlist.dat', 'ascii')
#events.gtis = [[events.time[0], events.time[1]], [events.time[150], events.time[151]]]
#gtis = np.asarray([[events.time[0], events.time[1]]])
gtis = np.asarray([[58000., 59000.], [59015., 59100.]])
#gtis = np.asarray([[58000., 59000.], [89015., 89100.]])
events.gtis = gtis
print(events.gtis)
print(events.time[:5])
#lc = Lightcurve.read('timeseries_0.dat', 'ascii')
#print(lc.time[:5], lc.counts[:5])

'''
plt.figure()
plt.plot(lc.time, lc.counts, color="blue")
plt.title("Light curve")
plt.xlabel(f"Time (s from {events.mjdref})")
plt.ylabel(f"Counts/bin")
plt.savefig('lightcurve.png')
plt.show()
'''

#'''
#folding for one period
f_profile = 'profile_mvm.png'
nbin = 32
test_periods = [10., 10.1, 10./3, 1.1]
fig, ax1 = plt.subplots()
#'''
#'''
# folding using events
for p in test_periods:
    ph, profile, profile_err = fold_events(events.time, 1/p, nbin=nbin, gti=gtis)
    ax1 = plot_profile(ph, profile, ax=ax1)
    fig.savefig(f_profile)


#'''

'''
# folding using lightcurve
for p in test_periods:
    ph, profile, profile_err = fold_events(lc.time, 1/p, nbin=nbin, weights=lc.counts)
    ax1 = plot_profile(ph, profile, ax=ax1)
    fig.savefig(f_profile)
'''

#'''
# We will search for pulsations over a range of frequencies around the known pulsation period.
nbin=10
df_min = 1/obs_length
oversampling=10
df = df_min / oversampling
frequencies = np.arange(1/period - 200 * df, 1/(period) + 200 * df, df)
#frequencies = np.arange(0.02, 1/(period/3) + 200 * df, df)
# entire observation
#frequencies = np.arange(0., obs_length, df)
#print(frequencies)
#'''
# filenames for plots
#profile = 'profile.png'
f_chi2 = 'chi2_mvm.png'
f_chi2_fit = 'chi2_fit_mvm.png'


# perform epoch folding on various test periods based on ###EVENTLIST###
freq, efstat = epoch_folding_search(events.time, frequencies, nbin=nbin, expocorr=True, gti=gtis)
###LIGHTCURVE###
#freq, efstat = epoch_folding_search(lc.time, frequencies, nbin=nbin, weights=lc.counts)
#print(lc.counts)
#'''

#'''
# search best candidate frequencies
ntrial = (frequencies[-1] - frequencies[0]) / df_min
ef_detlev = fold_detection_level(nbin, epsilon=0.001, ntrial=len(freq))
cand_freqs_ef, cand_stat_ef = search_best_peaks(freq, efstat, ef_detlev)
print(cand_freqs_ef)
# fit gaussian
fg=fit_gaussian(freq, efstat-(nbin-1),amplitude=max(efstat-(nbin-1)), 
                mean=cand_freqs_ef[0], stddev=1/(np.pi*obs_length))
# fit sinc
#fg=fit_sinc(freq, efstat-(nbin-1),amp=max(efstat-(nbin-1)), mean=cand_freqs_ef[0], 
#            obs_length=obs_length)
#print(fg.mean)
fit_f = fg.mean[0]
#'''

#'''
# ---- PLOTTING --------
plt.figure()
plt.plot(freq, efstat, label='EF statistics')
plt.axhline(nbin - 1, ls='--', lw=3, color='k', label='n - 1')
plt.axvline(1/period, lw=3, alpha=0.5, color='r', label='Correct frequency')
#plt.axvline(1/(period/3), lw=3, alpha=0.5, color='r', label='Correct frequency')
plt.xlabel('Frequency (Hz)')
plt.ylabel('EF Statistics')
_ = plt.legend()
plt.savefig(f_chi2)
#'''
#'''
# plot chi2 with fit
plt.figure(figsize=(15, 5))
plt.plot(freq, efstat-(nbin-1), label='EF statistics')
plt.plot(freq, fg(freq), label='Best fit')
plt.axvline(1/period, alpha=0.5, color='r', label=f'Correct frequency {1/period}')
#plt.axvline(1/(period/3), alpha=0.5, color='r', label=f'Correct frequency {1/(period*3)}')
plt.axvline(fg.mean[0], alpha=0.5, label=f'Fit frequency {fit_f}')
#plt.xlim( [fg.mean[0] * 0.95, fg.mean[0] * 1.05] )
plt.xlabel('Frequency (Hz)')
plt.ylabel('EF Statistics')
plt.legend()
plt.savefig(f_chi2_fit)
#'''
#'''

'''
#residuals
#plt.figure(figsize=(15, 5))
plt.plot(freq, efstat-(nbin-1)-fg(freq))
plt.xlabel('Frequency (Hz)')
_ = plt.ylabel('Residuals')

'''