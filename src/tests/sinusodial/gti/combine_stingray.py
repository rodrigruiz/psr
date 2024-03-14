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

#gtis = np.asarray([(58000., 59000.), (59015., 59100.)])
gtis = [(58000., 59000.)]
#gtis = [(58000., 59000.), (59015., 59100.)]
#gtis = np.asarray([[58000., 59000.], [59015., 59100.]])
#events.gtis = gtis
#print(events.gtis)
#print(events.time[:5])

ev1 = EventList().read('eventlist_0.dat', 'ascii')
ev1.gti = [(58000., 58010.)]
ev2 = EventList().read('eventlist_2.dat', 'ascii')
ev2.gti = [(59015., 59100.)]
events = ev1.join(ev2)
#events.gtis = gtis
print(ev1.gti, ev2.gti, events.gti)
events.write('combined_el2.dat', 'ascii')