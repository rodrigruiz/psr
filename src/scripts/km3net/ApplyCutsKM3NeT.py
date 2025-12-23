""" Apply selection cuts to the events from a KM3NeT event list, like source location or energy threshold. 

Usage: ApplyCutsKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR [--detector=<detector>] [--energy_min=<float>] [--energy_max=<float>] [--nhits_tr_min=<float>] [--nhits_sh_min=<float>] [--lik_tr_min=<float>] [--lik_sh_min=<float>] [--beta0_tr_max=<float>] [--zenith_min=<float>] [--trackscore_low_max=<float>] [--trackscore_high_min=<float>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files or file pattern
  -o --output_dir OUTPUT_DIR             Output directory
     --detector=<string>                 Detector location ('arca','orca','antares') [default: arca]
     --energy_min=<float>                Minimum energy in GeV [default: 0.]
     --energy_max=<float>                Maximum energy in GeV [default: 1e10]
     --nhits_tr_min=<float>              Minimum number of hits for track reconstructions [default: 0.]
     --nhits_sh_min=<float>              Minimum number of hits for shower reconstructions [default: 0.]
     --lik_tr_min=<float>                Minimum likelihood for track reconstructions [default: 0.]
     --lik_sh_min=<float>                Minimum likelihood for shower reconstructions, negated [default: 0.]
     --beta0_tr_max=<float>              Maximum beta0 value for track reconstructions [default: 1.0]
     --zenith_min=<float>                Minimum zenith value in degrees (Below horizon: >90deg) [default: 0.]
     --trackscore_low_max=<float>        Maximum trackscore to keep at the lower edge of the score distribution [default: 0.5]
     --trackscore_high_min=<float>       Minimum trackscore to keep at the higher end [default: 0.5]
"""
#python3 psr/src/scripts/CreateEventListKM3NeT.py -i '/home/hpc/capn/capn107h/software/hdf5TestOutput/*' -o eventlistTestOutput/ -s hdf5SourceFiles/Vela_X-1.h5

from docopt import docopt
import os, glob

from h5py import File, Group
from astropy.table import Table
from astropy.time import Time
from astropy.io.misc.hdf5 import read_table_hdf5, write_table_hdf5
import astropy.units as u
from astropy.coordinates import SkyCoord
import km3io.definitions as kd
import plens.EventList as EL
from km3astro.io import load_hdf5_tables

import numpy as np

from km3astro.coord import local_event
from km3astro import sources


def main():
    arguments = docopt(__doc__)
    #data = {key.replace("-", ""): arguments[key] for key in arguments}

    input_files = arguments['--input_files']
    output_dir = arguments['--output_dir']

    detector = str(arguments['--detector'])

    energy_min = float(arguments['--energy_min'])
    energy_max = float(arguments['--energy_max'])
    nhits_tr_min = float(arguments['--nhits_tr_min'])
    nhits_sh_min = float(arguments['--nhits_sh_min'])
    lik_tr_min = float(arguments['--lik_tr_min'])
    lik_sh_min = float(arguments['--lik_sh_min'])
    beta0_tr_max = float(arguments['--beta0_tr_max'])
    zenith_min = float(arguments['--zenith_min'])
    trackscore_low_max = float(arguments['--trackscore_low_max'])
    trackscore_high_min = float(arguments['--trackscore_high_min'])


    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    for file in input_files:
        with File(file, 'r') as h5_file:

            # Create Filename
            folder_path, file_name = os.path.split(file)
            file_name = os.path.splitext(file_name)[0]
            output_file = f"{arguments['--output_dir']}{file_name}_cutsapplied.hdf5"


            # Load File:            
            EventListNoCuts = EL.readEventList(h5_file)
        

            # Energy Cut:
            energy_mask = (EventListNoCuts['energy'].value >= energy_min) & (EventListNoCuts['energy'].value <= energy_max)
            print("Energy Mask: ", energy_mask)
            l_beforecut = len(EventListNoCuts)
            EventList = EventListNoCuts[energy_mask]
            l_aftercut = len(EventList)
            print(f"Cut on Energy applied: {energy_min} (E_min) and {energy_max} (E_max)")
            print(f"Events before Cut: {l_beforecut}")
            print(f"Events after Cut: {l_aftercut}")


            # Select Shower and Track Events:
            track_events = (EventList['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & (EventList['rec_stage'] == kd.reconstruction.JMUONBEGIN)
            if detector == 'arca':
                shower_events = (EventList['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & (EventList['rec_stage'] == kd.reconstruction.AASHOWERBEGIN)
                hits_mask_shower = shower_events & (EventList['AASHOWERFIT_NUMBER_OF_HITS'] >= nhits_sh_min)
            elif detector == 'orca':
                #shower_events = (EventList['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & (EventList['rec_stage'] == kd.reconstruction.AASHOWERBEGIN)
                shower_events = (EventList['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & (EventList['rec_stage'] == kd.reconstruction.JSHOWERBEGIN)
                hits_mask_shower = shower_events & (EventList['JGANDALF_NUMBER_OF_HITS'] >= nhits_sh_min)


            # Number of Hits mask
            hits_mask_track = track_events & (EventList['JGANDALF_NUMBER_OF_HITS'] >= nhits_tr_min)
            #hits_mask_shower = shower_events & (EventList['AASHOWERFIT_NUMBER_OF_HITS'] >= nhits_sh_min)
            hits_mask = hits_mask_shower | hits_mask_track
            print("Hits Mask: ", hits_mask)
            l_beforecut = len(EventList)
            EventList = EventList[hits_mask]
            l_aftercut = len(EventList)
            print(f"Cut on Hits applied with Thresholds: {nhits_tr_min} (track) and {nhits_sh_min} (shower)")
            print(f"Events before Cut: {l_beforecut}")
            print(f"Events after Cut: {l_aftercut}")

            # Select Shower and Track Events again:
            track_events = (EventList['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & (EventList['rec_stage'] == kd.reconstruction.JMUONBEGIN)
            if detector == 'arca':
                shower_events = (EventList['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & (EventList['rec_stage'] == kd.reconstruction.AASHOWERBEGIN)
            elif detector == 'orca':
                #shower_events = (EventList['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & (EventList['rec_stage'] == kd.reconstruction.AASHOWERBEGIN)
                shower_events = (EventList['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & (EventList['rec_stage'] == kd.reconstruction.JSHOWERBEGIN)

            # Likelihood mask
            lik_mask_track = track_events & (np.abs(EventList['likelihood']) >= lik_tr_min)
            lik_mask_shower = shower_events & (np.abs(EventList['likelihood']) >= lik_sh_min)
            lik_mask = lik_mask_track | lik_mask_shower
            print("Likelihood Mask: ", lik_mask)
            l_beforecut = len(EventList)
            EventList = EventList[lik_mask]
            l_aftercut = len(EventList)
            print(f"Cut on Likelihood applied with Thresholds: {lik_tr_min} (track) and {lik_sh_min} (shower)")
            print(f"Events before Cut: {l_beforecut}")
            print(f"Events after Cut: {l_aftercut}")

            # Zenith mask
            zenith_min_rad = zenith_min * np.pi/180
            zenith_mask = (EventList['zenith'] >= zenith_min_rad)
            print("Zenith Mask: ", zenith_mask)
            l_beforecut = len(EventList)
            EventList = EventList[zenith_mask]
            l_aftercut = len(EventList)
            print(f"Cut on Zenith applied with min zenith {zenith_min} degrees ({zenith_min_rad} rad)")
            print(f"Events before Cut: {l_beforecut}")
            print(f"Events after Cut: {l_aftercut}")

            # Select Shower and Track Events again:
            track_events = (EventList['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & (EventList['rec_stage'] == kd.reconstruction.JMUONBEGIN)
            if detector == 'arca':
                shower_events = (EventList['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & (EventList['rec_stage'] == kd.reconstruction.AASHOWERBEGIN)
            elif detector == 'orca':
                #shower_events = (EventList['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & (EventList['rec_stage'] == kd.reconstruction.AASHOWERBEGIN)
                shower_events = (EventList['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & (EventList['rec_stage'] == kd.reconstruction.JSHOWERBEGIN)

            # Beta0 Mask
            beta0_tr_mask = track_events & (EventList['JGANDALF_BETA0_RAD'] <= beta0_tr_max)
            beta0_mask = beta0_tr_mask | shower_events
            print("Beta0 Mask: ", beta0_mask)
            l_beforecut = len(EventList)
            EventList = EventList[beta0_mask]
            l_aftercut = len(EventList)
            print(f"Cut on Beta0 applied to track events with max beta0 {beta0_tr_max}")
            print(f"Events before Cut: {l_beforecut}")
            print(f"Events after Cut: {l_aftercut}")

            # Trackscore Mask
            trackscore_mask = (EventList['track_score'] <= trackscore_low_max) | (EventList['track_score'] >= trackscore_high_min)
            print("Trackscore Mask: ", trackscore_mask)
            l_beforecut = len(EventList)
            EventList = EventList[trackscore_mask]
            l_aftercut = len(EventList)
            print(f"Cut on Trackscore applied: Omitting events with track_score between {trackscore_low_max} and {trackscore_high_min}")
            print(f"Events before Cut: {l_beforecut}")
            print(f"Events after Cut: {l_aftercut}")


            # Store Meta Info 
            EventList.meta['energy_min'] = energy_min
            EventList.meta['energy_max'] = energy_max
            EventList.meta['nhits_tr_min'] = nhits_tr_min
            EventList.meta['nhits_sh_min'] = nhits_sh_min
            EventList.meta['lik_tr_min'] = lik_tr_min
            EventList.meta['lik_sh_min'] = lik_sh_min
            EventList.meta['beta0_tr_max'] = beta0_tr_max
            EventList.meta['zenith_min'] = zenith_min
            EventList.meta['trackscore_low_max'] = trackscore_low_max
            EventList.meta['trackscore_high_min'] = trackscore_high_min

            # Write HDF5 File
            write_table_hdf5(EventList, output_file, path='timeseries', overwrite =True, serialize_meta=True)
            print(f"All cuts applied. New file written to {output_file}")



if __name__ == "__main__":
    main()