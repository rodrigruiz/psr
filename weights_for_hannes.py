#!/usr/bin/env python

import os
import sys
import km3io as ki
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

#from km3astro.io import dir_to_spherical

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

use_orca = True

# ORCA
if use_orca:
    path = '/home/wecapstor3/capn/capn107h/orca/mcv9.2/'
    files = [path + f for f in os.listdir(path) if 'mc.gsg_neutrinos' in f and f.endswith('offline.v9.2.root')][:10]

# ARCA
else:
    path = '/home/wecapstor3/capn/capn107h/mc/'
    files = [path + f for f in os.listdir(path) if 'mcv8.1.gsg_' in f and '132' in f]


energies = np.zeros(shape=int(1e6))
cos_zeniths = np.zeros(shape=int(1e6))
weights = np.zeros(shape=int(1e6))
pids = np.empty(shape=int(1e6), dtype='S10')

livetime = 0

counter = 0
for file in files:

    f = ki.OfflineReader(file)

    mc_tracks = f.mc_tracks
    run_duration = f.header.DAQ.livetime
    energy = mc_tracks.E[:, 0]
    dir_z = mc_tracks.dir_z[:, 0]
    cos_zenith_true = -dir_z  # to be verified!
    #theta, phi = dir_to_spherical(mc_tracks[:,0])
    #azimuth, zenith = neutrino_to_source_direction(phi, theta)
    #azimuth_true = ki.tools.azimuth(mc_tracks[:,0])
    #print(azimuth)
    pdg = mc_tracks.pdgid[:, 0]

    w = f.w
    w2 = w[:, 1]

    # v9.2
    try:
        n_gen = 1 / w[:, 3]  # to be verified!
    # v8.1
    except IndexError:
        n_gen = f.header.genvol.numberOfEvents

    w2list = f.w2list
    is_cc = w2list[:, 10] == 2

    df = {'pdgid': pdg, 'is_cc': is_cc, 'n_gen': n_gen, 'w2': w2, 'energy': energy,
          'cos_zenith_true': cos_zenith_true}

    weight = wgts.get_weights_pid_output(df)
    pid = wgts.get_pid_class(pdg, is_cc)

    n_events = len(energy)
    energies[counter: counter + n_events] = energy
    cos_zeniths[counter: counter + n_events] = cos_zenith_true
    weights[counter: counter + n_events] = weight*run_duration
    pids[counter: counter + n_events] = pid
    counter += n_events

    livetime += f.header.DAQ.livetime

    f.close()

energies = energies[:counter]
cos_zeniths = cos_zeniths[:counter]
weights = weights[:counter]
pids = pids[:counter]
pid_set = list(set(pids))

if use_orca:
    bins = np.geomspace(1, 100, 51)
    label = 'orca'
else:
    bins = np.geomspace(100, 10000, 51)
    label = 'arca'

fig, ax = plt.subplots()

for p in pid_set:
    # if p != b'muon_cc':
    #     continue
    is_pid = pids == p
    # to see oscillations
    is_pid = np.logical_and(is_pid, cos_zeniths <= 0) # <=-0.2
    print(p, np.sum(is_pid), np.sum(weights[is_pid]))
    ax.hist(energies[is_pid], bins, weights=weights[is_pid], histtype='step', label=p)

ax.set_xscale('log')
ax.set_xlabel('energy [GeV]')
#ax.set_yscale('log')
#ax.set_ylabel('weighted event number [1/s]')
ax.set_ylabel('expected amount of events')

fig.legend()
fig.tight_layout()

fig.savefig('test_weights_%s_nevents.png' % label)

is_mu = pids == b'muon_cc'
z_bins = np.linspace(-1, 1, 41)

h, x, y = np.histogram2d(energies[is_mu], cos_zeniths[is_mu], bins=(bins, z_bins), weights=weights[is_mu])
fig, ax = plt.subplots()
xx, yy = np.meshgrid(x, y)
im = ax.pcolormesh(xx, yy, h.T) #, norm=mpl.colors.LogNorm())
plt.colorbar(im, ax=ax)
ax.set_xscale('log')
ax.set_xlabel('energy [GeV]')
ax.set_ylabel('cos_zenith')
fig.savefig('test_weights_2d_%s_nevents.png' % label)
