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
           
            if added_profile is None and err_buffer is None:
                added_profile = data[1][:]
                err_buffer = np.array(data[2][:])
            else:
                added_profile += data[1][:]
                err_buffer = np.column_stack((err_buffer, np.asarray(data[2][:]).T))
       
        added_profile_err = np.sqrt(np.sum((err_buffer**2), axis=1).tolist())

        return data[0][:], added_profile, added_profile_err

def read_hdf5(file):
    data = Table()
    data['phase_bins'] = file['ef_profile/phase_bins']
    data['profile'] = file['ef_profile/profile']
    data['profile_err'] = file['ef_profile/profile_err']
        
    return data
    
def add_profiles(input_dir, filepattern, frequency, outputdir):
    
    added_profile = None
    err_buffer = None
    
    input_files = glob.glob(input_dir + filepattern + str(frequency) + '.hdf5')
    #print(input_files)
    for file in input_files:
        
        with h5py.File(file) as f:
            data = read_hdf5(f)
           
        if added_profile is None and err_buffer is None:
            added_profile = data['profile']
            err_buffer = np.array(data['profile_err'])
        else:
            added_profile += data['profile']
            err_buffer = np.column_stack((err_buffer, np.asarray(data['profile_err']).T))
        
    added_profile_err = np.sqrt(np.sum((err_buffer**2), axis=1).tolist())
    
    added_data = { 'ef_profile/phase_bins' : data['phase_bins'], 'ef_profile/profile' : added_profile, 'ef_profile/profile_err' : added_profile_err }
    save_hdf5(outputdir+'Antares_added_profile_'+str(frequency)+'.hdf5', **added_data)
    
    return data['phase_bins'], added_profile, added_profile_err
