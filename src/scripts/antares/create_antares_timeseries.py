"""
Usage:
  create_antares_timeseries.py -f <input_files>... -e <efficiencies_file> [-h]


Options:
  -f <input_files>...    List of .hdf5 files with the ANTARES rates
  -e <efficiencies_file> Fits file with the summary table of ANTARES OM efficiencies as a function of time.
  -h                     Print help
"""
import antares.antares_hdf5 as ah5
from docopt import docopt
import h5py
from astropy.io import fits
from astropy.table import QTable
from astropy.time import Time, TimeDelta


def get_efficiency(t, eff_table):
    """
    Find the y value corresponding to the lower bound of the interval in the table for a given x.

    Parameters:
    - t: The input value for which to find the corresponding y.
    - eff_table: An Astropy QTable with columns 'unix_time' and 'mean_efficiency', sorted by 'unix_time'.

    Returns:
    - The 'mean_efficiency' value corresponding to the lower bound of the interval.
    """
    for row in range (1, len(eff_table)):
        low = 0
        high = len(eff_table)-1
        
        while low < high:

            middle = (low + high) //2

            if eff_table.iloc[middle]['unix_time'] < t.to_value('unix'):
                closest_index = middle
                low = middle + 1
            else:
                high = middle - 1

        return closest_index

def main():
    arguments = docopt(__doc__)

    efficiencies = QTable.read(arguments['-e'])
    efficiencies.sort('unix_time')
    efficiencies.add_index('unix_time')
    for file in arguments['-f']:
        run = h5py.File(file, 'r')
        unixtime = ah5.get_start_time(run)
        index = get_efficiency(unixtime, efficiencies)
        print(index)
        print(len(efficiencies))
        print(efficiencies[index]['year'], efficiencies[index]['month'],efficiencies[index]['day'])
        print(unixtime.to_value('datetime'))
        
if __name__ == '__main__':
    main()


