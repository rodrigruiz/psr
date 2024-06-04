""" Performs corrections on TimeSeries.

Usage: CorrectTimeSeries.py -i INPUT_FILE [--wildcard] -o OUTPUT_DIR -s SOURCE_SPECS_FILE

Options:
  -h --help                              Help
  -i --input_file INPUT_FILE             Input file
  -o --output_dir OUTPUT_DIR             Output file
  -s --source SOURCE_SPECS_FILE
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
from astropy.coordinates import EarthLocation, SkyCoord, Angle
from astropy.table import Table



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
    
    if data['wildcard']:
        file = glob.glob(data['input_file'][0])
    else:
        file = data['input_file']
    
    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    folder_path, file_name = os.path.split(data['input_file'])
    file_name = os.path.splitext(file_name)[0]

    output_file = data['output_dir'] + file_name + '_corrected.hdf5'
    print(output_file, os.path.exists(output_file))
    
        
    with h5py.File(data['source'], 'r') as f:
        # Read the source name
        source_name = f['source_name'][()]
        print("Source Name:", source_name)

        # Read the SkyCoord information
        ra = f['skycoord/ra'][()]
        dec = f['skycoord/dec'][()]
        skycoord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
        print("SkyCoord:", skycoord)

        # Read the orbital period
        Porb = f['orbital_period'][()]
        print("Orbital Period:", Porb, "days")
        
        Tpi2 = f['Tpi2'][()]
        axsini = f['axsini'][()]
        e = f['e'][()]
        omega = f['omega'][()]

    # Read TimeSeries
    with h5py.File(file) as input_file:
        
        TimeSeries = EL.readEventList(input_file)
        print("Raw TimeSeries:")
        print(TimeSeries)
        TimeSeries = EL.barycentric_correction(TimeSeries, skycoord)
        print("After Baryocentric Correction:")
        print(TimeSeries)
        time_corr = orbit_cor_deeter(TimeSeries.time.value, 
                                        (float(Porb)*u.d).to_value(u.s), 
                                        float(axsini), 
                                        float(e), 
                                        Angle(float(omega), u.deg).radian - np.pi/2, 
                                        Time(float(Tpi2) + float(Porb)/2, format='jd').unix
                                    )
        
        TimeSeries = Table([time_corr], names=['time'])
        print("After Binary Corrections:")
        print(TimeSeries)

    # Store Corrected TimeSeries again in HDF5-File
    with h5py.File(output_file, 'w') as output:  

        #print(TimeSeries)
        #EL.saveEventList(TimeSeries, output)
        TimeSeries.write(output, format='hdf5', overwrite=True, serialize_meta=True)
    
if __name__ == "__main__":
    main()
    