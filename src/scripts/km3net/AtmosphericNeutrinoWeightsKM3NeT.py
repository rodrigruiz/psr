""" Load KM3NeT root data files and convert them to astropy tables. 

Usage: AtmosphericNeutrinoWeightsKM3NeT.py -i INPUT_FILES... -o OUTPUT_DIR [--detector=<detector>]

Options:
  -h --help                              Help
  -i --input_files INPUT_FILES           Input files or file pattern
  -o --output_dir OUTPUT_DIR             Output directory
     --detector=<string>                 Detectorname ('orca' or 'arca') [default: 'arca']

"""

import sys
import os
import re
import numpy as np
import pandas as pd
from h5py import File
from astropy.table import Table
import matplotlib.pyplot as plt
import km3io as ki
import glob
from docopt import docopt

# Custom weights module path
sys.path.append('/home/hpc/capn/mppi083h/Code/work/')
import weights as wgts

def dir_to_spherical(tracks):
    """
    Convert directional coordinates (x, y, z) to spherical coordinates (theta, phi).

    Parameters
    ----------
    tracks : numpy.array
        Array containing track information, OfflineReader.tracks object

    Returns
    -------
    tuple
        Tuple containing spherical coordinates (theta, phi).
    """
    vec = np.array([tracks.dir_x, tracks.dir_y, tracks.dir_z])
    vec = np.transpose(vec)
    theta_val = ki.tools.theta(vec)
    phi_val = ki.tools.phi(vec)
    return theta_val, phi_val

def neutrino_to_source_direction(phi, theta, radian=True):
    """Flip the direction.

    Parameters
    ----------
    phi, theta: neutrino direction
    radian: bool [default=True]
        receive + return angles in radian? (if false, use degree)

    """
    phi = np.atleast_1d(phi).copy()
    theta = np.atleast_1d(theta).copy()
    if not radian:
        phi *= np.pi / 180
        theta *= np.pi / 180
    assert np.all(phi <= 2 * np.pi)
    assert np.all(theta <= np.pi)
    azimuth = (phi + np.pi) % (2 * np.pi)
    zenith = np.pi - theta
    if not radian:
        azimuth *= 180 / np.pi
        zenith *= 180 / np.pi
    return azimuth, zenith

def extract_run_id(filename: str) -> str:
    # Find all 8-digit numbers in the filename
    matches = re.findall(r'(?<!\d)\d{8}(?!\d)', filename)
    
    # Return the last occurrence if there are multiple
    if matches:
        return matches[-1]  # Last occurrence
    
    return None  # No 8-digit number found

def calculate_weights(input_files,detectorname='arca'):
    """
    Load event data, calculate weights, and save a new file with event_id and weight.
    """
    energies = np.zeros(shape=int(1e6))
    cos_zeniths = np.zeros(shape=int(1e6))
    new_weights = np.zeros(shape=int(1e6))
    event_ids = np.zeros(shape=int(1e6))
    run_ids = np.zeros(shape=int(1e6))
    pids = np.empty(shape=int(1e6), dtype='S10')

    counter = 0
    for input_file in input_files:
        # Read the event list

        folder_path, file_name = os.path.split(input_file)
        file_name = os.path.splitext(file_name)[0]

        run_id = extract_run_id(file_name)
        
        f = ki.OfflineReader(input_file)

        run_id_array = [run_id] * len(f)
        event_id = range(0, len(f))
        event_id_np = np.array(event_id)

        mc_tracks = f.mc_tracks
        run_duration = f.header.DAQ.livetime
        energy = mc_tracks.E[:, 0]
        dir_z = mc_tracks.dir_z[:, 0]
        cos_zenith_true = -dir_z  # to be verified!
        print("cos_zenith_true: ", cos_zenith_true)
        theta, phi = dir_to_spherical(mc_tracks[:,0])
        #print("theta: ",theta)
        #print("phi: ", phi)
        azimuth, zenith = neutrino_to_source_direction(phi, theta)
        #print("azimuth: ",azimuth)
        print("zenith: ", zenith)
        #print("cos zenith: ", np.cos(zenith))
        print("energy: ", energy)
        print("mean energy: ", np.mean(energy))

        pdg = mc_tracks.pdgid[:, 0]

        w = f.w
        w2 = w[:, 1]

        print("w2: ", w2)
        print("mean weight: ", np.mean(w2))
        # v9.2
        try:
            num_gen_events = 1 / w[:, 3]  # to be verified!
        # v8.1
        except IndexError:
            num_gen_events = f.header.genvol.numberOfEvents

        w2list = f.w2list
        is_cc = w2list[:, 10] == 2


        print("zenith: ", zenith)
        print("cos_zenith_true: ", cos_zenith_true)
        print("energy: ", energy)
        print("Mean energy: ", np.mean(energy))
        print("w2: ", w2)
        print("w2abs: ", w2)
        print("Mean weight: ", np.mean(w2))
        print("num_gen_events: ", num_gen_events)

        osc = False

        if detectorname == 'arca': 
            osc = False
        elif detectorname == 'orca': 
            osc = True
        else: print("Invalid detectorname!")

        df = {
            'pdgid': pdg,
            'is_cc': is_cc,
            'w2': w2,
            'n_gen': num_gen_events,
            'energy': energy,
            'cos_zenith_true': cos_zenith_true,
            'osc': osc
        }

        print(df)

        # Compute new weights
        new_weight = wgts.get_weights_pid_output(df)

        pid = wgts.get_pid_class(pdg, is_cc)


        print("new_weight: ", new_weight)
        print("Amount of events: ", len(new_weight))
        print("Mean of new_weights: ", np.mean(new_weight))
        print("Sum of new_weights: ", np.sum(new_weight))

        n_events = len(energy)
        energies[counter: counter + n_events] = energy
        cos_zeniths[counter: counter + n_events] = cos_zenith_true
        new_weights[counter: counter + n_events] = new_weight #*run_duration
        pids[counter: counter + n_events] = pid
        

        event_ids[counter: counter + n_events] = event_id_np
        run_ids[counter: counter + n_events] = run_id_array
        counter += n_events

        f.close()

    energies = energies[:counter]
    cos_zeniths = cos_zeniths[:counter]
    new_weights = new_weights[:counter]/counter
    event_ids = event_ids[:counter]
    run_ids = run_ids[:counter]
    pids = pids[:counter]
    pid_set = list(set(pids))
    if detectorname == 'arca':
        bins = np.geomspace(100, 1e8, 51)
    elif detectorname == 'orca':    
        bins = np.geomspace(1, 1000, 51)
    else: print(f"Invalid detectorname: '{detectorname}' - Must be arca or orca!")
    fig, ax = plt.subplots()

    for p in pid_set:
        is_pid = pids == p
        # to see oscillations
        #is_pid = np.logical_and(is_pid, cos_zenith_true <= 0) # <=-0.2
        print(p, np.sum(is_pid), np.sum(new_weights[is_pid]))
        ax.hist(energies[is_pid], bins, weights=new_weights[is_pid], histtype='step', label=p)

    ax.set_xscale('log')
    ax.set_xlabel('energy [GeV]')
    ax.set_yscale('log')
    #ax.set_ylabel('weighted event number [1/s]')
    ax.set_ylabel('expected amount of events')

    fig.legend()
    fig.tight_layout()

    fig.savefig('test_weights_%s_nevents.png' % detectorname)

    is_mu = pids == b'muon_cc'
    z_bins = np.linspace(-1, 1, 41)

    h, x, y = np.histogram2d(np.asarray(energies[is_mu]), np.asarray(cos_zeniths[is_mu]), bins=(bins, z_bins), weights=new_weights[is_mu])
    fig, ax = plt.subplots()
    xx, yy = np.meshgrid(x, y)
    im = ax.pcolormesh(xx, yy, h.T) #, norm=mpl.colors.LogNorm())
    plt.colorbar(im, ax=ax)
    ax.set_xscale('log')
    ax.set_xlabel('energy [GeV]')
    ax.set_ylabel('cos_zenith')
    fig.savefig('test_weights_2d_%s_nevents.png' % str(detectorname))

    #event_ids = range(0, counter)
    #event_ids_np = np.array(event_ids)

    # Create a DataFrame with event_id and new_weight
    # Event_id and Run_id for each event meed to be created!!
    output_df = pd.DataFrame({
        'run_id': run_ids,
        'event_id': event_ids,
        'new_weight': new_weights
    })

    run_ids = np.array(run_ids, dtype=int)
    event_ids = np.array(event_ids, dtype=int)
    new_weights = np.array(new_weights, dtype=float)
    detectorname = np.full_like(event_ids,detectorname,dtype='U10')

    table = Table(
        names=('run_id', 'event_id', 'new_weight','detector'),
        data=[run_ids, event_ids, new_weights, detectorname],
        dtype=('i8', 'i8', 'f8', 'U10') 
    )

    # Generate output file name
    base, ext = os.path.splitext(input_file)
    output_file = f"{base}_atm_neutrino_weights.hdf5"

    # Save the table to HDF5
    table.write(output_file, path='weights', format='hdf5', overwrite=True)

    print(f"Saved weights as Astropy table to {output_file}")


if __name__ == "__main__":
    arguments = docopt(__doc__)
    
    input_files = []
    for pattern in arguments['--input_files']:
        input_files.extend(glob.glob(pattern))
    input_files.sort()

    print(input_files)

    detectorname = str(arguments['--detector'])
    print("input_file_path: ", input_files)
    print("detectorname: ", detectorname)
    calculate_weights(input_files, detectorname)
