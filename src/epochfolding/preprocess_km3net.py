from h5py import File, Group
from astropy.table import Table
from astropy.time import Time
from astropy.io.misc.hdf5 import read_table_hdf5
import astropy.units as u
from astropy.coordinates import SkyCoord
import km3io.definitions as kd

from km3astro.io import load_hdf5_tables

import numpy as np

from km3astro.coord import local_event
from km3astro import sources

'''
def load_hdf5_tables(file_path):
    """
    Load an HDF5 file and assign tables to an object.
    The function expects a file that is generated with astropy_hdf5_from_root() and has the paths 'HEADER','ID','RECO_EVENTS' and optionally 'MC_EVENTS'

    Parameters
    ----------
    file_path : str
        Path to the HDF5 file.

    Returns
    -------
    object
        An object containing attributes header_table, id_table, reco_table, and optionally mc_table.
    """
    tables = type("Tables", (), {})()  # Creating an empty object

    try:
        with File(file_path, "r") as h5file:
            tables.header_table = read_table_hdf5(h5file, path="HEADER")
            tables.id_table = read_table_hdf5(h5file, path="ID")
            if "RECO/RECO_EVENTS" in h5file:
                print("RECO_EVENTS is present")
                tables.reco_table = read_table_hdf5(h5file, path="RECO/RECO_EVENTS")

            else:
                print("RECO_EVENTS is NOT present")
                # tables.reco_table = read_table_hdf5(h5file, path='RECO/RECO_EVENTS')
            if "MC/MC_EVENTS" in h5file:
                tables.mc_table = read_table_hdf5(h5file, path="MC/MC_EVENTS")
                print("MC_EVENTS is present")
            else:
                print("MC_EVENTS is NOT present")
                # tables.mc_table = read_table_hdf5(h5file, path='MC/MC_EVENTS')
                tables.mc_table = None
    except Exception as e:
        print(f"Error loading tables from {file_path}: {e}")

    return tables

    
'''

def eventlist_from_hdf5(filepath, reco_type = "jmuon", output = 'eventlist', format = 'ascii', detector_location = 'arca', source_location = None, energy_threshold = None, distance_deg = 10):
    tables = load_hdf5_tables(filepath)
    if reco_type == "mc":
        event_table = tables.mc_table
        times = tables.id_table['timeslice_utc_time']
        #print(len(times))
        #print(times)
    else:
        if reco_type == 'jmuon':
            reco_type_mask = (tables.reco_table['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & (tables.reco_table['rec_stage'] == kd.reconstruction.JMUONBEGIN)
        elif reco_type == 'aashower':
            reco_type_mask = (tables.reco_table['rec_type'] == kd.reconstruction.AANET_RECONSTRUCTION_TYPE) & (tables.reco_table['rec_stage'] == kd.reconstruction.AASHOWERBEGIN)
        elif reco_type == 'jshower':
            reco_type_mask = (tables.reco_table['rec_type'] == kd.reconstruction.JPP_RECONSTRUCTION_TYPE) & (tables.reco_table['rec_stage'] == kd.reconstruction.JSHOWERBEGIN)    
        elif reco_type == 'aashower':
            reco_type_mask = (tables.reco_table['rec_type'] == kd.reconstruction.DUSJ_RECONSTRUCTION_TYPE) & (tables.reco_table['rec_stage'] == kd.reconstruction.DUSJSHOWERBEGIN)
        
        event_table = tables.reco_table[reco_type_mask]
        times = event_table['tracktime_utc']
        #print(len(times))
        #print(times)
    
    if energy_threshold is not None:
        energy_mask = event_table['energy'] >= energy_threshold
        event_table = event_table[energy_mask]

    times = event_table['tracktime_utc'] if reco_type != "mc" else event_table['timeslice_utc_time']
    
    if source_location is not None:
        event_location = local_event(np.array(event_table['theta_detectorframe']),np.array(event_table['phi_detectorframe']),times,detector_location)
        separation = event_location.separation(source_location)
        print(f"separation: {separation}")
        location_mask = separation <= distance_deg * u.deg
        event_table = event_table[location_mask]

    times = event_table['tracktime_utc'] if reco_type != "mc" else event_table['timeslice_utc_time']
    print(f"reco_type: {reco_type}")
    print(f"energy_threshold: {energy_threshold}")

    event_list = Table([times],names=['time'])
    print(event_list)
    
    energy_threshold_str = str(energy_threshold) if energy_threshold is not None else "no"
    eventlist_filename = 'eventlistTestOutput/' + reco_type +  '-recotype_' + energy_threshold_str + '-energythreshold_' + output

    if format == 'hdf5':
        event_list.write(eventlist_filename + '.hdf5', format='hdf5', overwrite=True, serialize_meta=True)
    else:
        event_list.write(eventlist_filename + '.dat', format='ascii.ecsv', overwrite=True)


def load_source_hdf5(filepath):
    with File(filepath, 'r') as f:
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

    return source_name, skycoord, Porb, Tpi2, axsini, e, omega

hdf5_filepath = "hdf5TestOutput/mcv8.1.gsg_anue-CCHEDIS_1e2-1e8GeV.sirene.jterbr00013767.jchain.aashower.311.h5"
source_location = sources.VELA_X

#eventlist_from_hdf5(filepath=hdf5_filepath,reco_type='mc')
eventlist_from_hdf5(filepath=hdf5_filepath,reco_type='jmuon',format='hdf5',energy_threshold=None,source_location=source_location)
#eventlist_from_hdf5(filepath=hdf5_filepath,reco_type='aashower')