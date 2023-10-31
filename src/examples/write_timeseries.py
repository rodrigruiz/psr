from astropy.table import Table
from astropy.time import Time
from astropy.timeseries import TimeSeries
import numpy as np

# Generate synthetic time and flux data
time = Time(np.linspace(58000, 58500, 100000), format='mjd')
flux = np.random.normal(0, 0.1, 100000)

# Create a TimeSeries object with a single 'flux' column
ts = TimeSeries(flux, time_scale='utc')

# Write the TimeSeries to a FITS file
ts.write('timeseries.fits', format='fits', overwrite=True)
