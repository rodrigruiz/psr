import numpy as np
import pickle

def findGTIs(a, time, outputfile):
    """Finds 0. entries in data array and returns time intervals of good data.
    
    Parameters
    ----------
        a : np.array
            Data to test.
        
        time : np.array
            Times for the data.
            
    Returns
    -------
        gtis : np.array
            Good Time Intervals.
        
    """
    # Create an array that is 1 where a is 0, and pad each end with an extra 0.
    iszero = np.concatenate(([0], np.equal(a, 0).view(np.int8), [0]))
    absdiff = np.abs(np.diff(iszero))
    # Runs start and end where absdiff is 1.
    rate_zero = np.where(absdiff == 1)[0].reshape(-1, 2)
    #print(rate_zero)
    # correct 0 column entries for gtis
    for index, i in enumerate(rate_zero[:,0]):
        rate_zero[index][0] = i - 1
    #print(rate_zero)
    
    # 0th and last element
    indices = np.concatenate(([0], rate_zero.flatten(), [len(a)-1])).reshape(-1, 2)
    
    # check if first bin is bti
    if (indices[0] < 0.).any():
        indices = np.delete(indices, 0, 0)
    if (indices[-1] >= len(a)).any():
        indices = np.delete(indices, -1, 0)

    return list((map(tuple, time[indices])))

def saveGTIs(gtis, outputfile):
    """Save GTIs to pickle file.
    
    Parameters
    ----------
        gtis : np.array
            Good Time Intervals.
        
        outputfile : str
            Filename/Filpath to save GTIs to. Without extension.
            
    """
    
    with open(outputfile + '.pkl', 'wb') as f:
        pickle.dump(gtis, f)
        
    #with open('testpickle.pkl', 'rb') as f:
    #gtis_file = pickle.load(f)
    
    return

def loadGTIs(file):
    with open(file, 'rb') as f:
        gtis = pickle.load(f)

    return gtis