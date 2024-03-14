import numpy as np
from astropy.time import Time, TimeDelta
from astropy.timeseries import BinnedTimeSeries, TimeSeries
import h5py 
from stingray import Lightcurve
import pickle

def findGTIs(a, time):
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    rate_zero = np.where(absdiff == 1)[0].reshape(-1, 2)
    print(rate_zero)
    # correct 0 column entries for gtis
    for index, i in enumerate(rate_zero[:,0]):
        rate_zero[index][0] = i - 1
    print(rate_zero)
    
    # 0th and last element
    indices = np.concatenate(([0], rate_zero.flatten(), [len(a)-1])).reshape(-1, 2)
    
    # check if first bin is bti
    if (indices[0] < 0.).any():
        indices = np.delete(indices, 0, 0)
    if (indices[-1] >= len(a)).any():
        indices = np.delete(indices, -1, 0)
    print(indices)
    gtis = time[indices]
    
    return list((map(tuple, gtis)))



f = h5py.File('/home/wecapstor3/capn/mppi148h/total_rates/Antares_054420_total_rates_combined.hdf5', 'r')
ts = TimeSeries(time=Time(f['timeseries/time_bin_start'][()], format='unix') )
ts['rateOn'] = f['timeseries/rateOn'][()]


gtis = findGTIs(ts['rateOn'], ts.time.value)
#time = np.array([1., 2., 3., 4., 7.])
#test = np.array([1., 1., 1., 1., 1.])
#print(findGTIs(test, time))

#"""
print(type(gtis))
with open('testpickle.pkl', 'wb') as f:
        pickle.dump(gtis, f)
        
with open('testpickle.pkl', 'rb') as f:
    gtis_file = pickle.load(f)

#print(a)

lc = Lightcurve(ts.time.value, ts['rateOn'], dt=0.104858, skip_checks=True, gti=gtis_file)
#lc = Lightcurve(ts.time.value, ts['rateOn'], dt=0.104858, skip_checks=True, gti=gtis)
print(lc.gti[0][0])
#print(len(lc.time))
lc.apply_gtis()
#print(len(lc.time))
#"""