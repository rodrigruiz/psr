import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import TimeSeries
import astropy.units as u
from astropy.time import Time

# Generate a synthetic time series with noise and a periodic signal
np.random.seed(0)
n_points = 1000
time = np.linspace(0, 20, n_points)
periodic_signal = 2.0 * np.sin(2 * np.pi * time / 3.0)
noise = np.random.normal(0, 0.5, n_points)
flux = periodic_signal + noise

# Create a Time object and assign time data
time = Time(time, format='jd')  # Specify the time format

# Create a TimeSeries object and add the flux data
ts = TimeSeries(time=time)
ts['flux'] = flux

# Test different periods and create folded plots
test_periods = [2.0, 3.0, 4.0, 5.0]  # Replace with the periods you want to test

for period_value in test_periods:
    # Convert the period to an astropy Quantity with the correct units (e.g., days)
    period = period_value * u.day

    folded_data = ts.fold(period=period)

    plt.figure(figsize=(8, 6))
    plt.plot(folded_data['time'], folded_data['flux'], marker='o', linestyle='None')
    plt.xlabel('Folded Time (Phase)')
    plt.ylabel('Flux')
    plt.title(f'Folded Data (Period = {period_value} days)')
    plt.grid(True)

    plt.savefig(f'folded_data_period_{period_value}.png')
    plt.close()

# Show or save the plots for different test periods
plt.show()
