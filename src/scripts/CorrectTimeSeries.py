""" Performs corrections on TimeSeries.

Usage: CorrectTimeSeries.py -i INPUT_FILE -o OUTPUT_FILE [--rajd RIGHT_ASCENSION --decjd DECLINATION]

Options:
  -h --help                              Help
  -i --input_file INPUT_FILE             Input files
  -o --output_file OUTPUT_FILE           Output file
     --rajd RIGHT_ASCENSION              Right ascension (J2000) (degrees)
     --decjd DECLINATION                 Declination (J2000) (degrees)

"""
# python3 CorrectTimeSeries.py -i 'CombinedTimeSeries_Test.hdf5' -o 'CorrectedTimeSeries_Test.hdf5' --rajd '256.06166667' --decjd '-60.28166667'
# python3 CorrectTimeSeries.py -i "CombinedTimeSeries_MEff_MA.hdf5" -o "TimeSeries_J1704-6016.hdf5" --rajd "256.06166667" --decjd "-60.28166667"

from docopt import docopt

import h5py as h5py

from astropy.time import Time, TimeDelta
from astropy.timeseries import BinnedTimeSeries
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import EarthLocation


# PLENS Imports
import plens.TimeSeries as TS
import plens.antares_hdf5
import plens.antares_hdf5 as antares_hdf5




def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]

    # Read TimeSeries
    with h5py.File(data['input_file']) as input_file:
        TimeSeries, timeslice_duration = TS.readTimeSeries(input_file)  

    
    # Perform Barycentric Correction
    if data['rajd'] != None:
        print('Perform Barycentric Correction')
        skycoord = coord.SkyCoord(ra=data['rajd'], dec=data['decjd'], unit=(u.deg, u.deg))
        TimeSeries = TS.barycentric_correction(TimeSeries, timeslice_duration, skycoord)

    print('Barycentric Correction Successful')
    
    # Fill Gaps in TimeSeries
    TimeSeries = TS.fill_TimeSeries_gaps(TimeSeries, timeslice_duration)

    # Store Corrected TimeSeries again in HDF5-File
    with h5py.File(data['output_file'], 'w') as output_file:   
    
        TS.saveTimeSeries(TimeSeries, timeslice_duration, output_file)
        
    
if __name__ == "__main__":
    main()
    