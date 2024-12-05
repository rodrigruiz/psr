""" von Mises dist

Usage: PlotDistribution.py -o OUTPUT_DIR [--latex]

Options:
  -h --help                              Help
  <INPUT_FILES>...           Input files
  -o --output_dir OUTPUT_DIR             Output directory
     --nbin=<int>                        Number of bins in the folded profile [default: 32]
     --bins_hist=<int>                   Number of bins in the histogram [default: 200]
     --latex                             Whether to use latex formatting. [default: False]

"""
#  python3 antares_plot_chi2.py -i './folded/Antares_*_chi2_*.hdf5' -o './'
import os, glob
from docopt import docopt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from plens.plot_latex_size import set_size
from plens.PulseModel import MVMD

width = 418.25368
figsize = set_size(width, fraction=1, ratio='golden')

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]

def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

def main():
    arguments = docopt(__doc__)

    data = {}
    for key in arguments:
        data[key.replace("-", "")] = arguments[key]
    print(data)
    
    if data['latex']:
        plt.style.use('~/software/Psr/src/latex.mplstyle')
        
    
    output = data['output_dir'] + 'vonMises'
        
    fig, ax = plt.subplots(figsize=figsize)
    x = np.linspace(0, 10, 300)
    f=0.1
    kappa = [0.05, 5., 29., 57.]
    alpha = np.linspace(0.1, 1, len(kappa))

    
    #hex_list = ['#00008B', '#0000ff', '#ff0000', '#d62728']
    #hex_list = ['#0000ff', '#87cefa', '#d62728']
    #hex_list = ['#0000ff', '#ffffff', '#d62728']
    #hex_list = [ '#0000ff', '#b0c4de', '#87cefa', '#d3d3d3' ]
    #hex_list = [ '#0000ff', '#d3d3d3' ]
    hex_list = [ '#0000ff','#ffffff' ]
    #hex_list = [ '#0000ff', '#ff0000' ]
    #float_list = [0., 0.9, 1.]
    cm = get_continuous_cmap(hex_list)
    #colors = plt.cm.seismic(np.linspace(0.3,0.7,len(kappa)))
    #colors = plt.cm.Blues(np.linspace(0.2,1.,len(kappa)))
    #colors = plt.cm.coolwarm(np.linspace(0.,1.,len(kappa)))
    colors_cm = cm(np.linspace(0.,.9, len(kappa)))
    colors = ['#00008B', '#0000ff', colors_cm[2], colors_cm[3] ]
    print(colors)
    #for k, a in zip(kappa, alpha):
    #    plt.plot(x*f, MVMD(x, f, np.pi, k, 1.), alpha=a, color='b', label=r'$\kappa = {}$'.format(k))
    for k, c in zip(kappa, colors):
        plt.plot(x*f, MVMD(x, f, np.pi, k, 1.), color=c, label=r'$\kappa = {}$'.format(k))
    
    plt.legend()
    plt.title(r'modified von Mises distribution for various shape parameters $\kappa$')
    plt.xlabel(r'Phase $\varphi$')
    plt.ylabel('relative Intensity [a.u.]')
    fig.savefig(output + '.pdf',  bbox_inches='tight')
    #plt.show()
    plt.close()
        
if __name__ == "__main__":
    main()