""" Performs corrections on TimeSeries.

Usage: FFT_Analysis.py -i INPUT_FILE -o OUTPUT_FILE --jname JNAME [-f FREQUENCY] [-H HORDER] [-l] [-b]

Options:
  -h --help                              Help
  -i --input_file INPUT_FILE             Input files
  -o --output_file OUTPUT_FILE           Output file
  -j --jname JNAME                       Pulsar name based on J2000 coordinates
  -f --frequency FREQUENCY               Pulsar frequency
  -H --horder HORDER                     Highes Harmonic Order (possible values [2,4,8])
  -l --logscale                          Add logscale plot of Natural Spectrum
  -b --blind                             Blind the data
  
"""
# python3 FFT_Analysis.py -i 'CorrectedTimeSeries_Test.hdf5' -o 'Test.txt'

from docopt import docopt

import h5py as h5py

import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import norm, chi2

from astropy.timeseries import BinnedTimeSeries
from astropy.table import vstack

import plens.TimeSeries as TS
import plens.RedNoiseFilter as RNF
import plens.HarmonicSum as HS

import re






def sigma(n):
    # Get fraction within n standard deviations for a normal distribution
    return norm.cdf(n) - norm.cdf(-1*n)


def power_above_expectation(spectrum, df=2, H=1):
    # Power value above the expected number of sample points is <1
    
    pow_a_e = HS.chi2_HS(df, H).isf(1/len(spectrum))
    return pow_a_e


def create_subplot(x, y, ax, color, title, semilogx=False, semilogy=False, sigma_axis=False, df=2, H=1, frequency=None):
    
    # Plot Data
    ax.plot(x, y, color=color)
    
    # Set Title and Labels
    ax.set_title(title)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Power [a.u.]')
    ax.grid()
    
    # Logscale
    if semilogx == True:
        ax.semilogx()
#    else:
#        ax.set_xlim([0,5])
    if semilogy == True:
        ax.semilogy()

    # Add second y-axis on the right in sigma values
    
    n = np.arange(11)# Integers nums of sigmas
    
    if sigma_axis == True:
        
        # alpha: significance level (Need to divide by number of points in the spectrum, due to look-elsewhere effect
        N = len(y)
        alpha = ( 1 - sigma(n) ) / N

        ticks = list(HS.chi2_HS(df,H).isf(alpha)) + [power_above_expectation(y,df,H)]
        labels = [str(i)+r'$\sigma$' for i in n] + [r'$E$']
        
        ax2 = ax.secondary_yaxis('right',functions=(lambda x: x, lambda x: x))
        ax2.set_yticks(ticks=ticks,labels=labels)
        
        ax.axhline(power_above_expectation(y,df,H), color='tab:grey', linestyle='dashed')
    
    if frequency != None:
        ax.axvline(x=float(frequency), color='tab:red', alpha=0.8, zorder=1)
    
    return


def find_candidate(freq, spectrum, rate_key, sensitivity, df=2, H=1):
    # Find candidates in the spectrum
    
    # Calculate p-value
    p_val = HS.chi2_HS(df=df, H=H).sf(spectrum)
    
    candidate_info = list()
    
    for candidate, frequency, power, p in zip(p_val <= sensitivity, freq, spectrum, p_val):
        if candidate:
            print('Candidate found: %s \tFrequency = %f Hz,\tPower = %f,\tp-val = %f,\tH = %d' %(rate_key, frequency, power, p, H))
            candidate_info.append(np.array([rate_key, frequency, power, p, H]))
    
    return candidate_info




def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    
    if data['horder']:
        data['horder'] = int(data['horder'])
    
    # Read TimeSeries
    with h5py.File(data['input_file']) as input_file:
        TimeSeries, timeslice_duration = TS.readTimeSeries(input_file)

    print('The TimeSeries has a length of %d (2^%.1f) data points' %(len(TimeSeries), np.log2(len(TimeSeries))))    
        
    # Cut Dataset for Master Theses
    from astropy.time import Time
    TimeSeries = TimeSeries[TimeSeries['time_bin_start'] > Time('2010-11-01')]
    
    print('The TimeSeries has a length of %d (2^%.1f) data points after 2010-11-01' %(len(TimeSeries), np.log2(len(TimeSeries))))
    
    ##################################################    
    # Shuffle TimeSeries to blind data
    ##################################################    
    if data['blind']:
        
        # Random Number Generator
        rng = np.random.default_rng(seed=0)
        
        for key in TimeSeries.keys():
            rng.shuffle(TimeSeries[key].value)    
        
        
    ##################################################    
    # Expand TimeSeries to the next length
    ##################################################
    N_TS_total = 2**27
    N_add =  N_TS_total - len(TimeSeries)
    
    TimeSeries = vstack([TimeSeries,
                         BinnedTimeSeries(time_bin_start=TimeSeries['time_bin_start'][-1] + timeslice_duration,
                                          time_bin_size=timeslice_duration, 
                                          data={key: np.mean(TimeSeries[key]) * 
                                                np.ones(N_add, dtype=TimeSeries[key].dtype)
                                                for key in TimeSeries.keys() 
                                                if key not in ['time_bin_start', 'time_bin_size']}
                                         )
                        ])
    
    print('TimeSeries was extended by %d (2^%.1f) points to 2^%d points' %(N_add, np.log2(N_add), np.log2(N_TS_total)))
    
    
    ##################################################
    # Calculate FFT
    ##################################################
    rfftfreq  = np.fft.rfftfreq(len(TimeSeries), d=timeslice_duration.sec)
    rfftpower = {key : np.abs(np.fft.rfft(TimeSeries[key].value, norm='ortho'))**2
                for key in TimeSeries.keys()
                if key not in ['time_bin_start', 'time_bin_size']}
    
    # Degrees of freedom of the underlying chi-squared distribution
    df = 2

    
    ##################################################
    # Apply RedNoiseFilter
    ##################################################
    N_total = 2**26 # Total Number of points in the spectrum (without 0'th frequency bin) (N_total = N_TS_total/2)
    N_1 = 2**19 # Number of points in the first segment size
    n_1 = 16 # Number segments of size N_1
    N_2 = 2**20 # Number of points in the second segment size
    n_2 = 56 # Number of segments of size N_2
    # Note: N_total = N_1*n_1 + N_2*n_2

    segments = np.concatenate([
        np.linspace(N_1, N_1*n_1, n_1, endpoint=True, dtype=int),
        np.linspace(N_1*n_1 + N_2, N_total, n_2, endpoint=True, dtype=int)
    ])[:-1]
    # The values in 'segments' are the indices of the (inclusive) lower borders of each segment,
    # the border at 0 Hz is not needed ,therefore the first np.linspace starts at N_1
    # the border at the Nyquist frequency is not needed, therefore the slicing [:-1] of the concatenated array
    
    RNFpower = {key : RNF.RedNoiseFilter(rfftpower[key][1:], segments,  exclude_max=True)
                 for key in rfftpower.keys()
                }
    
    
    ##################################################
    # Apply HarmonicSummation
    ##################################################
    if data['horder'] != None:
    
        H_HS = 2**np.linspace(1, np.log2(data['horder']), int(np.log2(data['horder'])), dtype=int) # Highes Harmonic Orders to sums
        HSpower = {key : {h : HS.HarmonicSum(RNFpower[key][1:], h) 
                          for h in H_HS} 
                   for key in rfftpower.keys()}
        # Dictionary containing a Dictionary
        # First index is 'rateOff'/'rateOn', second index is H the highes summed Harmonic Order

    
    ##################################################
    # Create Plots of Spectra
    ##################################################    
    num_rows = 2 + data['logscale']
    if data['horder'] != None:
        num_rows += int(np.log2(data['horder']))
    # 2 plots for Natural Spectrum and RNF Spectrum,
    # if data['logscale'] == True, then data['logscale'] == 1 (as True == 1)
    # add one plot for horder=2, two for horder=4, and three for horder=8

    #fig, ax = plt.subplots(num_rows,2, figsize=(3.5*(3/2)*2, 3.5*num_rows)) # Factor (3/2) as the hight to widht ratio, the factor 2 due to two plots per row
    
    # inch to mm
    mm = 1/2.54 * 0.1
    fig, ax = plt.subplots(num_rows,2, figsize=(210*mm*5/4, 297*mm*5/4 * num_rows/6 * 0.95)) 

    fig.suptitle(data['jname'])

    rates = ['rateOff', 'rateOn']
    
    # Define Colors
    color = {'rateOff' : 'tab:blue',
             'rateOn'  : 'tab:orange'
             }
    # Defina Axis ID
    ax_id = {'rateOff' : 0,
             'rateOn'  : 1
            }
    
    # Create first rateOff plot, then rateOn
    for rate_type in rates:
        row = 0 # run variable to iterate over the rows
    
        # Natural Power Spectrum lin-log-scale
        create_subplot(rfftfreq[1:], rfftpower[rate_type][1:],
                       ax[row, ax_id[rate_type]],
                       color=color[rate_type], title=('%s\nNatural Spectrum' %rate_type), semilogy=True, frequency=data['frequency'])
        
        ax[row, ax_id[rate_type]].set_xlim([-0.1,5])
        
     
        # Natural Power Spectrum lin-log-scale
        if data['logscale']:
            row += 1
            
            create_subplot(rfftfreq[1:], rfftpower[rate_type][1:],
                           ax[row, ax_id[rate_type]],
                           color=color[rate_type], title=('Natural Spectrum'), semilogx=True, semilogy=True, frequency=data['frequency'])
            
            #ax[row, ax_id[rate_type]].set_xlim([-0.1,5])
                    
        # RedNoiseFilter Spectrum
        row +=1
        create_subplot(rfftfreq[1+N_1:], RNFpower[rate_type][N_1:], 
                       ax[row, ax_id[rate_type]],
                       color=color[rate_type], title='Red Noise Filter', sigma_axis=True, frequency=data['frequency'])
        
        if data['horder'] != None:
            for h in H_HS:
                row += 1
                # HSpower H = h
                create_subplot(rfftfreq[1+N_1:len(HSpower[rate_type][h])+1], HSpower[rate_type][h][N_1:],
                               ax[row, ax_id[rate_type]],
                               color=color[rate_type], title=(r'Harmonic Sum $H=%d$' %h), sigma_axis=True, H=h, frequency=data['frequency'])        
            
                ax[row, ax_id[rate_type]].set_xlim([-0.1,5])
    
    plt.tight_layout()
    plot_file_name = re.sub('\.\w*', '.pdf', data['output_file']) 
    plt.savefig(plot_file_name)
    
    
    
    ##################################################
    # Find Candidates
    ##################################################
    candidate_info = list()
    
    for key in RNFpower.keys():
        candidate_info += find_candidate(rfftfreq[1+N_1:], RNFpower[key][N_1:],
                                         key, 1/len(RNFpower[key]), df=df, H=1
                                        )
    if data['horder'] != None:
        for key in HSpower.keys():
            for h in HSpower[key].keys():
                candidate_info += find_candidate(rfftfreq[1+N_1:len(HSpower[key][h])+1], HSpower[key][h][N_1:],
                                                 key, 1/len(HSpower[key][h]), df=df, H=h)

    np.savetxt(data['output_file'], candidate_info, 
               header='Periodic Candidates\nRateType Frequency [Hz] Power [a.u.] p_val [] H []', fmt='%s')

    
if __name__ == "__main__":
    main()