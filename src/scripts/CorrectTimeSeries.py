""" Performs corrections on TimeSeries.

Usage: CorrectTimeSeries.py -i INPUT_FILES -o OUTPUT_DIR --filepattern=<filepattern> --correction=<correction> [--rajd RIGHT_ASCENSION --decjd DECLINATION] [--Porb=<float> --axsini=<float> --e=<float> --omega=<float> --Tpi2=<float>] (--timeseries | --eventlist)

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output file
     --correction=<correction>           Which correction to perform on the TimeSeries. Options: 'bary', 'fill', 'bary+bin', 'bary+fill' [default: bary+bin]
                                         If 'bary' or 'bary+bin' or 'bary+fill' is selected RA and Dec should also be specified.
                                         If 'bary+bin' is selected, additionally, Porb, axsini, e, omega and Tpi2 have to be specified. 
                                         If known these parameters are reported in https://gammaray.nsstc.nasa.gov/gbm/science/pulsars.html
                                         Please give the values of the parameters in the units given on the website.
     --rajd RIGHT_ASCENSION              Right ascension (J2000) (degrees)
     --decjd DECLINATION                 Declination (J2000) (degrees)
     
"""
#python3 CorrectTimeSeries.py -i '../stingray/antares/data/Antares_051870_total_rates_*[!combined].hdf5' -o '../stingray/antares/data/' --correction bary --rajd 02h43m40.4252869512s --decjd +61d26m03.757456824s
# python3 CorrectTimeSeries.py -i 'CombinedTimeSeries_Test.hdf5' -o 'CorrectedTimeSeries_Test.hdf5' --rajd '256.06166667' --decjd '-60.28166667'
# python3 CorrectTimeSeries.py -i "CombinedTimeSeries_MEff_MA.hdf5" -o "TimeSeries_J1704-6016.hdf5" --rajd "256.06166667" --decjd "-60.28166667"

from docopt import docopt
import os, glob
import h5py as h5py
import re
import numpy as np
from astropy.time import Time, TimeDelta
from astropy.timeseries import BinnedTimeSeries
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import EarthLocation
from astropy.coordinates import Angle


# PLENS Imports
import plens.TimeSeries as TS
import plens.EventList as EL
from plens.TimeSeries import orbit_cor_deeter
import plens.antares_hdf5
import plens.antares_hdf5 as antares_hdf5




def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    if data['rajd'] != None and data['decjd'] != None:
        skycoord = coord.SkyCoord(ra=data['rajd'], dec=data['decjd'], unit=(u.deg, u.deg))
    else:
        raise ValueError(
            "please specify RA and Dec of the object!")
    
    input_files = glob.glob(data['input_files'])
    input_files.sort()
    #print(input_files)
    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])
        
    # Construct TimeSeries for each run
    for file in input_files:
        
        #split = re.split('Antares_(\d*)_total_rates_(\d*).hdf5', file)
        split = re.split(data['filepattern'], file)
        run_number = split[1]
        split_number = split[2]
            
        print('Processing Run Nr.: ' + str(run_number) + ', split: ' + str(split_number))#, end='\r')
        
        
            
        # Read TimeSeries
        with h5py.File(file) as input_file:
            
            if data['timeseries']:
                TimeSeries, timeslice_duration = TS.readTimeSeries(input_file)  
                arg = '_total_rates_'
            elif data['eventlist']:
                TimeSeries = EL.readEventList(input_file)
                arg = '_eventlist_'
            
            #print(TimeSeries)
            output_file = data['output_dir'] + 'Antares_' + run_number + arg + split_number + '_corrected.hdf5'
            #print(output_file, os.path.exists(output_file))
            if os.path.exists(output_file):
                continue
            # Perform Barycentric Correction
            if data['correction'] == 'bary':
           
                #print('Perform Barycentric Correction')
                TimeSeries = TS.barycentric_correction(TimeSeries, skycoord)
                #print('Barycentric Correction Successful')
            
        
            elif data['correction'] == 'bary+bin':
                #print(data['eventlist'])
                if data['timeseries']:
                    print(TimeSeries.time_bin_start[0])
                    # barycentric correction
                    TimeSeries = TS.barycentric_correction(TimeSeries, skycoord)
                    print(TimeSeries.time_bin_start[0])
                    # binary correction
                    #TimeSeries['time_bin_start'].format = 'mjd'
                    # use tdb format; convert Tpi2 to tdb or unix
                    time_corr = orbit_cor_deeter(TimeSeries.time_bin_start.value, 
                                             (float(data['Porb'])*u.d).to_value(u.s), 
                                             float(data['axsini']), 
                                             float(data['e']), 
                                             Angle(float(data['omega']), u.deg).radian - np.pi/2, 
                                             Time(float(data['Tpi2']) + float(data['Porb'])/2, format='jd').unix
                                            )
                    TimeSeries['time_bin_start'] = time_corr
                    print(TimeSeries.time_bin_start[0])
                elif data['eventlist']:
                    #print(TimeSeries.time[0])
                    # barycentric correction
                    TimeSeries = EL.barycentric_correction(TimeSeries, skycoord)
                    #print(TimeSeries.time[0])
                    # binary correction
                    #TimeSeries['time_bin_start'].format = 'mjd'
                    # use tdb format; convert Tpi2 to tdb or unix
                    time_corr = orbit_cor_deeter(TimeSeries.time.value, 
                                             (float(data['Porb'])*u.d).to_value(u.s), 
                                             float(data['axsini']), 
                                             float(data['e']), 
                                             Angle(float(data['omega']), u.deg).radian - np.pi/2, 
                                             Time(float(data['Tpi2']) + float(data['Porb'])/2, format='jd').unix
                                            )
                    TimeSeries['time'] = time_corr
                    #print(TimeSeries.time[0])
                
            elif data['correction'] == 'fill':
                # Fill Gaps in TimeSeries
                TimeSeries = TS.fill_TimeSeries_gaps(TimeSeries, timeslice_duration)
        
            elif data['correction'] == 'bary+fill':
        
                #print('Perform Barycentric Correction')
                TimeSeries = TS.barycentric_correction(TimeSeries, skycoord)
                #print('Barycentric Correction Successful')
        
                # Fill Gaps in TimeSeries
                TimeSeries = TS.fill_TimeSeries_gaps(TimeSeries, timeslice_duration)
    
            else:
                raise ValueError(
                    "correction was not specified or the correction was not recognised.")
    
        # Store Corrected TimeSeries again in HDF5-File
        with h5py.File(output_file, 'w') as output:  
            if data['timeseries']:
                TS.saveTimeSeries(TimeSeries, timeslice_duration, output)
            elif data['eventlist']:
                #print(TimeSeries)
                #EL.saveEventList(TimeSeries, output)
                TimeSeries.write(output, format='hdf5', overwrite=True, serialize_meta=True)
    
if __name__ == "__main__":
    main()
    