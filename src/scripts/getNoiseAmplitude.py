# Get Noise Amplitude of the TimeSeries
import h5py as h5py
import numpy as np
import matplotlib.pyplot as plt
import plens.TimeSeries as TS
from datetime import datetime
from scipy.stats import chi2

input_file = "/home/saturn/capn/mppi107h/Data/antares_analysis/CombinedTimeSeries.hdf5"
#Test
#input_file = "/home/saturn/capn/mppi107h/Data/antares_analysis/TimeSeries_Test.hdf5"


with h5py.File(input_file) as file:

    TimeSeries, timeslice_duration = TS.readTimeSeries(file)
    
##################################################
# Due to mistake in the HDF5-File exstraction script, rateOn values need to be divided by 2
TimeSeries['rateOn'] = TimeSeries['rateOn']/2


##################################################    
    

mean   = "Mean \t %.2f \t %.2f"   %(np.mean(  TimeSeries['rateOff'].value)/1e3, np.mean(  TimeSeries['rateOn'].value)/1e3)
median = "Median \t %.2f \t %.2f " %(np.median(TimeSeries['rateOff'].value)/1e3, np.median(TimeSeries['rateOn'].value)/1e3)
std    = "Std \t %.2f \t %.2f"    %(np.std(   TimeSeries['rateOff'].value)/1e3, np.std(   TimeSeries['rateOn'].value)/1e3)
stats =  mean + "\n" + median + "\n" + std
print('in kHz\t rateOff \t rateOn')
print(stats)


##################################################
# Create Plot
##################################################

rate_types = ['rateOff', 'rateOn']
color = {'rateOff' : 'tab:blue','rateOn'  : 'tab:orange'}

mm = 1/25.4
fig = plt.figure(figsize=(145*mm, 145*mm * 4/6))
gs = fig.add_gridspec(1, 2,  width_ratios=(4, 1),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.1, hspace=0.05)

ax = fig.add_subplot(gs[0, 0])
ax_histy = fig.add_subplot(gs[0, 1], sharey=ax)
ax_histy.tick_params(axis="y", labelleft=False)



#ax.plot( (TimeSeries['time_bin_start'].unix - TimeSeries['time_bin_start'][0].unix)/(60*60*24), TimeSeries['rateOff'].value/1e3, 
#        marker='o', linestyle='None', markersize=1,  markeredgecolor='None', color=color['rateOff'], label='rateOff')
#ax.plot( (TimeSeries['time_bin_start'].unix - TimeSeries['time_bin_start'][0].unix)/(60*60*24), TimeSeries['rateOn'].value/1e3, 
#        marker='o', linestyle='None', markersize=1,  markeredgecolor='None', color=color['rateOn'], label='rateOn')

times = TimeSeries['time_bin_start'].to_value(format='datetime64')
bins = 100

ax.plot(times, TimeSeries['rateOff'].value/1e3,
        marker='o', linestyle='None', markersize=1,  markeredgecolor='None', color=color['rateOff'], label='rateOff')
ax.plot(times, TimeSeries['rateOn'].value/1e3,
        marker='o', linestyle='None', markersize=1,  markeredgecolor='None', color=color['rateOn'], label='rateOn')

    
ax_histy.hist(TimeSeries['rateOn'].value/1e3, bins=bins, density=True, histtype='stepfilled', orientation='horizontal',
              facecolor=color['rateOn'], edgecolor=color['rateOn'], alpha=0.75, range=(0, TimeSeries['rateOn'].max().value/1e3))

ax_histy.hist(TimeSeries['rateOff'].value/1e3, bins=bins, density=True, histtype='stepfilled', orientation='horizontal',
              facecolor=color['rateOff'], edgecolor=color['rateOff'], alpha=0.75, range=(0, TimeSeries['rateOff'].max().value/1e3))



ax.tick_params(axis='x', labelrotation=15)

ax.set_xlim([np.datetime64("2010-09-01").astype(datetime),np.datetime64("2011-02-14").astype(datetime)])

ax.set_ylabel('Rate [kHz]')
ax_histy.set_xlabel('Probability')
ax.legend(loc='upper right')

plt.savefig("Figures/TimeSeriesRate.png", dpi=300, bbox_inches = "tight")


##################################################
