#!/usr/bin/env python3
from astropy.timeseries import TimeSeries
from astropy.time import Time
import numpy as np

times = Time(np.linspace(58000,58500,100000), format='mjd')
periodic_signal = 0.5 * np.sin(2 * np.pi * times.jd / 10.0)
noise = np.random.normal(0, 0.1, 100000)
flux = periodic_signal + noise

ts = TimeSeries(time=times)
ts['flux']=flux
ts.write('timeseries.dat',format='ascii')
