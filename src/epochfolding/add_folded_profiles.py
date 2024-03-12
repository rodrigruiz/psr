#!/usr/bin/env python3
"""
Adding folded profiles.

Usage:
  add_folded_profiles.py --directory=<directory> --frequencies=<frequencies>
  add_folded_profiles.py (-h | --help)

Options:
  -h --help                   Show this help message and exit.
"""

import os, glob
import numpy as np
import h5py
from docopt import docopt
from astropy.table import Table
from astropy.io import ascii
from epochfolding.stingray_epochfolding import save_to_ascii, save_hdf5

def read_table(filename):
    data = ascii.read(filename, format='ecsv', fast_reader=True)
    return data
    
    
def add_folded_profiles(directory, frequencies):
    rootdir = os.path.normpath('./folded_events/')
    dir_list = [d for d in glob.glob(rootdir + '**/*') if os.path.isdir(d)]
    #print(dir_list)
    for p in frequencies:
        added_profile = None
        err_buffer = None
        for diry in dir_list:
            if os.path.isfile(diry + f'/folded_profile_{p}.dat'):
                file = diry + f'/folded_profile_{p}.dat'
            else:
                print(f"{diry} + /folded_profile_{p}.dat does not exist.")
                break 
            data = read_table(file)
            print(data[0])
            #print(data[0][:])
            if added_profile is None and err_buffer is None:
                added_profile = data[1][:]
                err_buffer = np.array(data[2][:])
            else:
                added_profile += data[1][:]
                err_buffer = np.column_stack((err_buffer, np.asarray(data[2][:]).T))
        #print(err_buffer.shape, len(added_profile))
        #print(err_buffer[0])
        added_profile_err = np.sqrt(np.sum((err_buffer**2), axis=1).tolist())
        #print(added_profile_err[0], len(added_profile_err))
        #save_to_ascii(data[0][:], added_profile, added_profile_err, p, outputdir= './added_profiles/', output='added_profiles_')
        return data[0][:], added_profile, added_profile_err

def read_hdf5(file):
    data = Table()
    data['phase_bins'] = file['ef_profile/phase_bins']
    data['profile'] = file['ef_profile/profile']
    data['profile_err'] = file['ef_profile/profile_err']
        
    return data
    
def add_profiles(input_dir, filepattern, frequency, outputdir):
    #rootdir = os.path.normpath('./folded_events/')
    #dir_list = [d for d in glob.glob(rootdir + '**/*') if os.path.isdir(d)]
    #for p in frequencies:
    print(frequency)
    added_profile = None
    err_buffer = None
    #data = []
    #for diry in dir_list:
    input_files = glob.glob(input_dir + filepattern + str(frequency) + '.hdf5')
    print(input_files)
    for file in input_files:
        #if os.path.isfile(diry + '/folded_profile_{0:.4f}.dat'.format(frequency)):
        #    file = diry + '/folded_profile_{0:.4f}.dat'.format(frequency)
        #else:
         #   print("{0}/folded_profile_{1:.4f}.dat does not exist.".format(diry, frequency))
        #    break 
        #data = read_table(file)
        with h5py.File(file) as f:
            data = read_hdf5(f)
            print(data[0])
        if added_profile is None and err_buffer is None:
            added_profile = data['profile']
            err_buffer = np.array(data['profile_err'])
            #print(err_buffer.shape)
            #print(added_profile)
            #print(err_buffer)
        else:
            added_profile += data['profile']
            err_buffer = np.column_stack((err_buffer, np.asarray(data['profile_err']).T))
            print(err_buffer.shape)
        #print(added_profile)
        #print(err_buffer)
        #print(err_buffer.shape, len(added_profile))
        #print(err_buffer[0])
        
    added_profile_err = np.sqrt(np.sum((err_buffer**2), axis=1).tolist())
    
    print(added_profile[0], added_profile_err[0], len(added_profile_err))
    #save_to_ascii(data[0][:], added_profile, added_profile_err, frequency, outputdir= './added_profiles/', output='added_profiles_')
    
    added_data = { 'ef_profile/phase_bins' : data['phase_bins'], 'ef_profile/profile' : added_profile, 'ef_profile/profile_err' : added_profile_err }
    save_hdf5(outputdir+'Antares_added_profile_'+str(frequency)+'.hdf5', **added_data)
    
    return data['phase_bins'], added_profile, added_profile_err
        
if __name__ == '__main__':
    arguments = docopt(__doc__)
    directory = arguments['--directory']
    if arguments['--frequencies'][0] != '[':
        frequencies = float(arguments['--frequencies'])
    else: 
        print(arguments['--frequencies'].strip('[]').split(','))
        frequencies = list( map(float, arguments['--frequencies'].strip('[]').split(',')) )
        print(list(frequencies))
    add_folded_profiles(directory, frequencies)
