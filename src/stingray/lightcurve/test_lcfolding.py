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

def sinusoid(times, frequency, baseline, amplitude, phase):
    return baseline + amplitude * np.sin(2 * np.pi * (frequency * times + phase))


period = 1.203501
mean_countrate = 50
pulsed_fraction = 0.2
bin_time = 0.01
obs_length = 3000

t = np.arange(0, obs_length, bin_time)
#print(t[0])

# The continuous light curve
counts = sinusoid(t, 1 / period, mean_countrate, 
                  0.5 * mean_countrate * pulsed_fraction, 0) * bin_time
lc = Lightcurve(t, counts, gti=[[-bin_time / 2, obs_length + bin_time / 2]],
                dt=bin_time)

plt.figure()
plt.plot(t[:100], counts[:100], color='blue' )
plt.xlabel('Time [a.u.]')
plt.ylabel('Counts [a.u.]')
plt.title('Sinusoid')
plt.grid(True)
#plt.legend()
plt.savefig( 'sine_stingray.png' )
plt.show()

'''
nbin = 32

ph, profile, profile_err = fold_events(events.time, 1/period, nbin=nbin)
_ = plot_profile(ph, profile)

ph, profile, profile_err = fold_events(events.time, 1/1.1, nbin=nbin)
_ = plot_profile(ph, profile)
'''