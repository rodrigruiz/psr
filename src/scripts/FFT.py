""" FFT

Usage: FFT.py -o OUTPUT_DIR [--latex]

Options:
  -h --help                              Help
  <INPUT_FILES>...           Input files
  -o --output_dir OUTPUT_DIR             Output directory
     --latex                             Whether to use latex formatting. [default: False]

"""
from docopt import docopt
import numpy as np
from plens.PulseModel import MVMD, sinusoid
from scipy.fft import fft, rfft
from scipy.fft import fftfreq, rfftfreq
from plens.plot_latex_size import set_size
import matplotlib.pyplot as plt
width = 418.25368
figsize = set_size(width, fraction=1, ratio='golden')

def plot_pulsetrain(times, counts, output='pulsetrain', color='blue', label=None, title='PulseTrain', fraction=1):
    #plt.figure(figsize=(5,3))
    plt.figure(figsize=set_size(width,fraction=fraction))
    plt.plot(times, counts, label=label, color=color )
    plt.xlabel('Time $t$ [s]')
    #plt.xlabel('Time [MJD]')
    plt.ylabel('Counts [a.u.]')
    plt.title(title)
    plt.grid(True)
    if label is not None:
        plt.legend()
    plt.savefig( output + '.pdf' )
    plt.show()
    
def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    print(data)
    
    if data['latex']:
        plt.style.use('~/software/Psr/src/latex.mplstyle')
        
    output = data['output_dir'] + 'FFT'
    
    length = 1000
    sampling_rate = 10
    times = np.linspace(0.,length, length*sampling_rate, endpoint=False)
    f, P_0, kappa, a = 0.2, 0, 5., 1.
    flux = MVMD(times, f, P_0, kappa, a)
    #f, baseline, amplitude, phase = 0.1, 0., 1., 0.
    #flux = sinusoid(times, f, baseline, amplitude, phase)
    plot_pulsetrain(times[:1000], flux[:1000], output='pulsetrain_mvm_zoom', fraction=0.495, title=r'$f_\mathrm{{MVMD}}(t; f, \kappa)$')
    plot_pulsetrain(times, flux, output='pulsetrain_mvm', fraction=0.495, title=r'$f_\mathrm{{MVMD}}(t; f, \kappa)$')
    noise = np.random.poisson(max(flux)*3,length*sampling_rate)
    signal = flux + noise
    plot_pulsetrain(times[:1000], signal[:1000], output='pulsetrain+noise_mvm_zoom', fraction=0.495, title=r'$f_\mathrm{{MVMD}}(t; f, \kappa)$')
    plot_pulsetrain(times, signal, output='pulsetrain+noise_mvm', fraction=0.495, title=r'$f_\mathrm{{MVMD}}(t; f, \kappa)$')
    
    #fourier = fft(signal)
    N = len(signal)
    
    #print(rfftfreq(N, d=1/sampling_rate), 2*np.abs(rfft(signal))/N)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    plt.plot(rfftfreq(N, d=1/sampling_rate)[1:], 2*np.abs(rfft(signal)[1:])/N, color='b')
    #plt.axvline(f, color='#d62728', ls='--', label='$f_\mathrm{signal} = 0.1\,\mathrm{Hz}$')
    #plt.legend()
    plt.title(r'FFT of $f_\mathrm{{MVMD}}(t; f = {} \,\mathrm{{Hz}}, \kappa = {})$'.format(f, kappa))
    #plt.title(r'FFT of a sinusoidal pulse train $ f(t; f = {} \,\mathrm{{Hz}})$'.format(f) + '$ = \sin{(2 \pi f t)} $')
    plt.xlim(0., 1.)
    plt.xlabel(r'Frequency $f$ [Hz]')
    plt.ylabel('Power spectral density')
    fig.savefig(output + '.pdf',  bbox_inches='tight')
    #plt.show()
    plt.close()
        


if __name__ == "__main__":
    main()