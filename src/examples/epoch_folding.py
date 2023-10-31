import numpy as np
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

# Generate synthetic time series data with a known periodic signal
np.random.seed(0)
t = np.linspace(0, 100, 1000)
# Inject a sinusoidal signal with a period of 10 seconds (0.1 Hz)
signal_period = 10
signal_amplitude = 5
signal = signal_amplitude * np.sin(2 * np.pi * (1 / signal_period) * t)
noise = np.random.normal(0, 1, len(t))
data = signal + noise

# Candidate frequencies to test (reciprocal of periods)
candidate_frequencies = 1 / np.linspace(5, 15, 1000)

# Initialize arrays to store periodogram results
periods = []
power = []

# Perform epoch folding and calculate the periodogram for each candidate frequency
for freq in candidate_frequencies:
    folded_time = t % (1 / freq)
    ls = LombScargle(folded_time, data)
    frequency, power_density = ls.autopower()

    # Find the frequency corresponding to the highest peak
    best_frequency = frequency[np.argmax(power_density)]

    # Convert frequency to period and store the results
    best_period = 1 / best_frequency
    periods.append(best_period)
    power.append(np.max(power_density))

# Plot the periodogram
plt.figure(figsize=(10, 4))
plt.plot(periods, power)
plt.xlabel('Period (units of time)')
plt.ylabel('Power')
plt.title('Periodogram')
plt.grid(True)

# Find the period corresponding to the highest peak in the periodogram
best_period = periods[np.argmax(power)]
print(f"Detected Period: {best_period:.2f} units of time")

plt.show()
