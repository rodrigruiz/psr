import os, glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
import plens.TimeSeries as TS

with h5py.File('table_links.h5',mode='w') as h5fw:
    link_cnt = 0 
    for h5name in glob.glob('file*.h5'):
        link_cnt += 1
        h5fw['link'+str(link_cnt)] = h5py.ExternalLink(h5name,'/')