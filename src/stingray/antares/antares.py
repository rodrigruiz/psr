import h5py
import numpy as np
import matplotlib.pyplot as plt
import plens
#import sys
#sys.path.insert(0, '../../plens/')

# do analysis on 0-ending runs only!!!
with h5py.File('/home/wecapstor3/capn/mppi19/ANTARES/data/rates/hdf5/Antares_051870_total_rates.hdf5', 'r') as f:
    print("Keys: %s" % f.keys())
    #print(f['lcmID_dict'])
    #print(f['header'].keys())
    print(f['total_rates'].attrs)