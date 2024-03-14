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
from docopt import docopt
from astropy.table import Table
from astropy.io import ascii
from stingray_epochfolding import save_to_ascii

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
        print(added_profile_err[0], len(added_profile_err))
        save_to_ascii(data[0][:], added_profile, added_profile_err, p, outputdir= './added_profiles/', output='added_profiles_')
        return data[0][:], added_profile, added_profile_err
    
def add_profiles(directory, frequency):
    rootdir = os.path.normpath('./folded_events/')
    dir_list = [d for d in glob.glob(rootdir + '**/*') if os.path.isdir(d)]
    #print(dir_list)
    #for p in frequencies:
    added_profile = None
    err_buffer = None
    for diry in dir_list:
        if os.path.isfile(diry + '/folded_profile_{0:.4f}.dat'.format(frequency)):
            file = diry + '/folded_profile_{0:.4f}.dat'.format(frequency)
        else:
            print("{0}/folded_profile_{1:.4f}.dat does not exist.".format(diry, frequency))
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
    print(added_profile_err[0], len(added_profile_err))
    save_to_ascii(data[0][:], added_profile, added_profile_err, frequency, outputdir= './added_profiles/', output='added_profiles_')
    
    return data[0][:], added_profile, added_profile_err
        
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
