from astropy.io import fits
import glob
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
plt.style.use('~/software/Psr/src/latex.mplstyle')
import numpy as np
#plt.style.use(['science', 'nature'])
#plt.rcParams['figure.dpi'] = 250
np.random.seed(20170615)
from tatpulsar.utils.functions import met2mjd, mjd2met
from plens.TimeSeries import orbit_cor_deeter#, orbit_cor_kepler
#from tatpulsar.pulse.binary import orbit_cor_kepler
from tatpulsar.pulse.barycor.barycor import barycor
from astropy.time import Time, TimeDelta
import astropy.units as u
""" to run in shell
# Download the Example Data of IXPE on Cen X-3
!wget -q -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np -R 'index*' -erobots=off --retr-symlinks \
https://heasarc.gsfc.nasa.gov/FTP/ixpe/data/obs/01//01006501/event_l2/
!wget -nH --no-check-certificate --cut-dirs=5 -r -l0 -c -N -np -R 'index*' -erobots=off --retr-symlinks \
https://heasarc.gsfc.nasa.gov/FTP/ixpe/data/obs/01/01006501/hk/ixpe01006501_all_orb_v07.fits.gz

!ls 01006501/*
"""

## READ Data
file = glob.glob("ixpe01006501_det*.fits")
orbit = glob.glob("ixpe01006501_all_orb_v07.fits")
time = np.array([])

for evt in file:
    print(evt)
    hdulist = fits.open(evt)
    time = np.append(time, hdulist[1].data.field("TIME"))

##EF
from stingray.pulse.search import epoch_folding_search
period = 4.79660
nbin = 64
df_min = 1e-2
df = 1e-6
frequencies = np.arange(1/period - df_min, 1/period + df_min, df)

freq, efstat = epoch_folding_search(met2mjd(time, telescope='ixpe'), frequencies, nbin=nbin)

fig, ax = plt.subplots()
#plt.rcParams['figure.dpi'] = 250
#plt.figure()
plt.plot(freq, efstat)
plt.axhline(nbin -1, ls='--', lw=3, color='k', label='n-1')
plt.axvline(1/period, lw=3, alpha=0.5, color='r', label='True frequency')
plt.axvline(freq[np.argmax(efstat)], lw=3, alpha=0.5, color='k', label='Observed frequency')
plt.xlabel("Freq (Hz)")
plt.ylabel("EF stat")
plt.legend()
fig.savefig('EF_wo_corrections.png', bbox_inches='tight')



ra=170.313301 
dec=-60.6233

# our barycentric correction
#from astropy import time, coordinates as coord, units as u

""" doesnt work; lihgt_travel_time works just for a location on earth
skycoord = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
times_original = time.Time(np.array(met2mjd(time, telescope='ixpe')), format='mjd')
delta = times_original.light_travel_time(skycoord=skycoord, kind='barycentric', location=????)
tdb_mjd = times_original + delta
tdb_met = mjd2met(tdb_mjd, 'ixpe')
"""
tdb_mjd = barycor(met2mjd(time, telescope='ixpe'),
                  ra=ra, dec=dec,
                  orbit=orbit[0],
                  accelerate=True)
tdb_met = mjd2met(tdb_mjd, 'ixpe')

##EF
freq, efstat = epoch_folding_search(tdb_met, frequencies, nbin=nbin)

fig, ax = plt.subplots()
#plt.rcParams['figure.dpi'] = 250
#plt.figure()
plt.plot(freq, efstat)
plt.axhline(nbin -1, ls='--', lw=3, color='k', label='n-1')
plt.axvline(1/period, lw=3, alpha=0.5, color='r', label='True frequency')
plt.axvline(freq[np.argmax(efstat)], lw=3, alpha=0.5, color='k', label='Observed frequency')
plt.xlabel("Freq (Hz)")
plt.ylabel("EF stat")
plt.legend()
fig.savefig('EF_bary_correction.png', bbox_inches='tight')

###BINARY CORRECTION###
#from tatpulsar.pulse.binary import orbit_cor_deeter

#-----------------------------------
Porb = 2.0869953 #day
axsini = 39.653 #light-sec
e = 0.0000
omega = 0.00 - np.pi/2
T_pi2 = 2455073.68504 - 2400000.5 #mjd 
x_mid = 298
y_mid = 306

Tnod = T_pi2 + Porb/2 # mjd
print(Tnod)
tnod = mjd2met(Tnod, telescope='ixpe') #met
print(tnod)
#time_mjd = met2mjd(time, telescope='ixpe')
print(tdb_mjd[0], type(tdb_mjd), type(tdb_mjd[0]))
print((tdb_mjd[0] - Tnod)*86400, tdb_met[0] - tnod)
#time = Time(time, format='jd').mjd
#print(time)
#print(time.jd)
#t = [t.value for t in time]
#tdb_bicor = orbit_cor_deeter(tdb_met, Porb*86400, axsini, e, omega, tnod)
#tdb_bicor = orbit_cor_deeter(tdb_mjd*86400, Porb*86400, axsini, e, omega, Tnod*86400)
tdb_bicor = orbit_cor_deeter(Time(tdb_mjd, format='mjd').unix, (Porb*u.d).to_value(u.s), axsini, e, omega, Time(Tnod, format='mjd').unix)

## EF search
freq, efstat = epoch_folding_search(tdb_bicor, frequencies, nbin=nbin)

print( f'actual frequency: {1/period}, observed: {freq[np.argmax(efstat)]}; difference: {np.abs(1/period - freq[np.argmax(efstat)])}')
# Plotting
fig, ax = plt.subplots()
plt.plot(freq, efstat)
plt.axhline(nbin -1, ls='--', lw=3, color='k', label='n-1')
plt.axvline(1/period, lw=3, alpha=0.5, color='black', label='True frequency')
plt.axvline(freq[np.argmax(efstat)], lw=3, alpha=0.5, color='r', label='Observed frequency')
plt.xlabel("Freq (Hz)")
plt.ylabel("EF stat")
plt.legend()
fig.savefig('EF_all_corrected_deeter.png', bbox_inches='tight')

"""
##Kepler
T_periastron = T_pi2 + Porb/2 # mjd
t_periastron = mjd2met(Tnod, telescope='ixpe') #met

tdb_bicor_kepler = orbit_cor_kepler(met2mjd(tdb_met, 'ixpe'),
                                    Tw=T_periastron,
                                    ecc=e,
                                    Porb=Porb*86400,
                                    omega=omega,
                                    axsini=axsini)
tdb_bicor_kepler = mjd2met(tdb_bicor_kepler, 'ixpe')

freq, efstat = epoch_folding_search(tdb_bicor_kepler, frequencies, nbin=nbin)

# Plotting
fig, ax = plt.subplots()
plt.plot(freq, efstat)
plt.axhline(nbin -1, ls='--', lw=3, color='k', label='n-1')
plt.axvline(1/period, lw=3, alpha=0.5, color='r', label='True frequency')
plt.axvline(freq[np.argmax(efstat)], lw=3, alpha=0.5, color='k', label='Observed frequency')
plt.xlabel("Freq (Hz)")
plt.ylabel("EF stat")
plt.legend()
fig.savefig('EF_all_corrected_kepler.png', bbox_inches='tight')
"""