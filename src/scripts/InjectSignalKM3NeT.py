"""Inject an artificial signal.
Usage: InjectSignal.py -i INPUT_FILES... -o OUTPUT_DIR [--ratio=<float>] [--pulseshape=<pulseshape>] [--df=<float>] [--frequency=<float>] [--baseline=<float>] [--a=<float>] [--phi=<float>] [--kappa=<float>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files
  -o --output_dir OUTPUT_DIR             Output file
     --ratio=<float>                     Ratio of injected signal to original count number. [default: 0.3]
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
from plens.PulseModel import MVMD, sinusoid

# PLENS Imports
import plens.EventList as EL
import plens.antares_hdf5
import plens.antares_hdf5 as antares_hdf5

from stingray import EventList, Lightcurve

def injectSignalRedistribute( time, ratio, bin_time, pulseshape, frequency, baseline, a, phi, kappa=None ):

    """Injects a signal with a MVM pulseshape into an existing time sequence (eventlist) by first deleting random counts to ensure an unchanged total count number.
    
    Parameters
    ----------
        time : np.array

        ratio : float
        
        bin_time : float
        
        frequency : float
            Frequency of the pulse train.
            
        baseline : float
            Offset along the y-axis.
            
        a : float
            Amplitude of the Pulse. Equates to the area of one pulse.
            
        phi : float
            Phase offset of the pulse train.
            
        kappa : float
            Shape parameter giving the width of the function.
        
    Returns
    -------
        np.array
        New times with injected pulsetrain.
        
    """
    
    if pulseshape == 'mvm':
        counts = MVMD(time, frequency, phi, kappa, a, baseline=baseline)
    elif pulseshape == 'sine':
        counts = sinusoid(time, frequency, baseline, a, phi)

    print(f"time: {time}")
    print(f"counts: {counts}")
    print(f"len(counts): {len(counts)}")

    # Create a light curve with the desired signal
    lc = Lightcurve(time, counts, dt=bin_time, skip_checks=True)

    lc.plot()
    
    # Simulate event times from the light curve
    ev = EventList()
    ev.simulate_times(lc)
    #print(f"len(ev.time): {len(ev.time)}")

    new_event_times = ev.time
    
        # Determine the desired number of new events based on the ratio
    total_events = len(time)
    num_new_events = int(ratio * total_events)
    num_original_events_to_keep = total_events - num_new_events

    print(f"len(time): {len(time)}")

    print(f"num_new_events: {num_new_events}")
    print(f"num_original_events_to_keep: {num_original_events_to_keep}")

    
    if num_original_events_to_keep < 0:
        raise ValueError("Number of original events to keep is negative. "
                         "Ensure that ratio is between 0 and 1.")
    
    # Randomly select events to remove from the new events if necessary
    if len(new_event_times) > num_new_events:
        indices_to_keep_new = np.random.choice(len(new_event_times), num_new_events, replace=False)
        new_event_times = new_event_times[indices_to_keep_new]
    
    # Randomly select events to keep from the original times
    indices_to_keep_original = np.random.choice(len(time), num_original_events_to_keep, replace=False)
    remaining_original_events = time[indices_to_keep_original]
    
    #print(f"len(new_event_times): {len(new_event_times)}")
    #print(f"len(remaining_original_events): {len(remaining_original_events)}")
    
    # Combine the remaining original events with the new events
    combined_events = np.sort(np.concatenate((remaining_original_events, new_event_times)))

    return combined_events

def main():
    #print(stingray.__version__)
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    input_files = []
    for pattern in data['input_files']:
        input_files.extend(glob.glob(pattern))
    input_files.sort()

    if not os.path.exists(data['output_dir']):
        os.makedirs(data['output_dir'])

    for file in input_files:    
        # Fetching filename for usage in output filename 
        folder_path, file_name = os.path.split(file)
        file_name = os.path.splitext(file_name)[0]
        #print(file_name)

        output_file = data['output_dir'] + file_name + '_'+ str(data['frequency']) + 'Hz_' + str(data['a']) + '-signal'
        #print(output_file, os.path.exists(output_file))

        # Read TimeSeries
        with h5py.File(file) as input_file:
            
            EventList = EL.readEventList(input_file)
            
            if data['pulseshape'] == 'mvm':
                EventListNew = TimeSeries(time=Time(injectSignalRedistribute(EventList['time'].value, 
                                                        float(data['ratio']),
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
                EventListNew = TimeSeries(time=Time(injectSignalRedistribute(EventList['time'].value, 
                                                        float(data['ratio']),
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
                
        if os.path.exists(output_file):
            print("File already existed. Deleting File...")
            os.remove(output_file)  # Remove the file if it already exists
        with h5py.File(output_file, 'w') as output:  

            EventListNew.write(output, format='hdf5', overwrite=True, serialize_meta=True)

if __name__ == "__main__":
    main()
