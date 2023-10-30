import numpy as np
from astropy.time import Time
from astropy.table import Table
from astropy.timeseries import TimeSeries
from astropy.io import fits

# Generate synthetic time and flux data with noise and a periodic signal
n_points = 100
time = Time(np.linspace(59000, 59100, n_points), format='mjd')

# Create a periodic signal (sine wave)
periodic_signal = 0.5 * np.sin(2 * np.pi * time.jd / 10.0)

# Add random noise to the periodic signal
noise = np.random.normal(0, 0.1, n_points)

# Combine the periodic signal and noise to create the flux
flux = periodic_signal + noise

# Create a table with time and flux columns
data_table = Table([time, flux], names=['time', 'flux'])

# Create a TimeSeries object from the table
timeseries = TimeSeries(data_table)

# Define the filename for the FITS file
fits_filename = 'timeseries_data_with_signal.fits'

# Save the TimeSeries to a FITS file
timeseries.write(fits_filename, format='fits', overwrite=True)

print(f'TimeSeries data with noise and a periodic signal saved to {fits_filename}')
