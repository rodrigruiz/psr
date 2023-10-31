from astropy.timeseries import TimeSeries
import matplotlib.pyplot as plt
ts = TimeSeries.read('timeseries.dat', format='ascii', time_column = 'time', time_format = 'mjd')
plt.plot(ts.time.mjd, ts['flux'], 'k.', markersize = 1)
plt.xlabel('Julian Date')
plt.ylabel('flux')
plt.show()
