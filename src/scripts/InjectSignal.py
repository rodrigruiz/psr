"""Inject an artificial signal.
Usage: InjectSignal.py -i INPUT_FILES -o OUTPUT_DIR --filepattern=<filepattern> [--pulseshape=<pulseshape>] [--df=<float>] [--frequency=<float>] [--baseline=<float>] [--a=<float>] [--phi=<float>] [--kappa=<float>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output file
     --pulseshape=<pulseshape>           Shape of the injected signal (sine or mvm). [default: mvm]
                                         if 'sine': 'df', 'frequency', 'baseline', 'a', 'phi' should be set
     --df=<float>                        Time resolution of the signal. [default: 0.1]
     --frequency=<float>                 Frequency of the signal. [default: 1]
     --baseline=<float>                  Offset on the y-axis. [default: 0.]
     --a=<float>                         Amplitude of the signal. [default: 1.]
     --phi=<float>                       Phase of the signal. [default: 0.]
     --kappa=<float>                     Shape parameter of the MVMD. [default: 5.]
"""

from docopt import docopt
import os, glob
import h5py as h5py
import re
import numpy as np
from astropy.time import Time, TimeDelta
from astropy.timeseries import TimeSeries
import astropy.units as u

# PLENS Imports
import plens.EventList as EL
import plens.antares_hdf5
import plens.antares_hdf5 as antares_hdf5

def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    input_files = glob.glob(data['input_files'])
    input_files.sort()

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])
        
    for file in input_files:
        
        #split = re.split('Antares_(\d*)_total_rates_(\d*).hdf5', file)
        split = re.split(data['filepattern'], file)
        run_number = split[1]
        split_number = split[2]
            
        print('Processing Run Nr.: ' + str(run_number) + ', split: ' + str(split_number))#, end='\r')
        
        output_file = data['output_dir'] + 'Antares_' + run_number + '_eventlist_' + split_number + '_signal'
        
        if os.path.exists(output_file):
            continue   
            
        # Read TimeSeries
        with h5py.File(file) as input_file:
            
            EventList = EL.readEventList(input_file)
            
            if data['pulseshape'] == 'mvm':
                EventListNew = TimeSeries(time=Time(EL.injectSignal(EventList['time'].value, 
                                                     float(data['df']), 
                                                     data['pulseshape'],
                                                     float(data['frequency']), 
                                                     float(data['baseline']), 
                                                     float(data['a']), 
                                                     float(data['phi']), 
                                                     float(data['kappa'])),
                                                format='unix'
                                                   )
                                         )
                output_file += '_mvm.hdf5'
                    
            elif data['pulseshape'] == 'sine':
                EventListNew = TimeSeries(time=Time(EL.injectSignal(EventList['time'].value, 
                                                     float(data['df']), 
                                                     data['pulseshape'],
                                                     float(data['frequency']), 
                                                     float(data['baseline']), 
                                                     float(data['a']), 
                                                     float(data['phi'])),
                                                format='unix'
                                                   )
                                         )
                output_file += '_sine.hdf5'
                
        with h5py.File(output_file, 'w') as output:  

            EventListNew.write(output, format='hdf5', overwrite=True, serialize_meta=True)

if __name__ == "__main__":
    main()
    