import numpy as np
from sys import exit
import os, sys
from copy import deepcopy
import pickle

def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out

def mergeDictionaries(genie_dict, nugen_dict, nugen_ccfactor, nugen_ncfactor, weight_keys, break_energy):
    # First correct NUGEN with a factor
    for weight_key in ['weight'] + weight_keys:
        nugen_dict['CC'][weight_key] /= nugen_ccfactor
        nugen_dict['NC'][weight_key] /= nugen_ncfactor

    skip_keys = ['GENIE_W']

    # Then merge the keys:
    new_dict = {}
    for interaction_key in genie_dict.keys():
        new_dict[interaction_key] = {}
        for event_info_key in genie_dict[interaction_key].keys():
            if event_info_key in skip_keys:
                continue
            genie_bool = genie_dict[interaction_key]['energy'] < break_energy
            nugen_bool = nugen_dict[interaction_key]['energy'] > break_energy
            new_dict[interaction_key][event_info_key]=\
                np.concatenate((genie_dict[interaction_key][event_info_key][genie_bool], 
                                nugen_dict[interaction_key][event_info_key][nugen_bool]),
                               axis = 0)
    
    return new_dict

def checkDictParams(default_parameters, user_parameters, warning_origin = None):
    wrong_keys = False
    user_keys = user_parameters.keys()
    for one_key in user_keys:
        if not default_parameters.has_key(one_key) and one_key!="sin_sq_th23":
            print warning_origin, ': Setting ', one_key, ' was given but is not used. Did you write it right?'
            wrong_keys = True

    if wrong_keys:
        print 'Check the input parameters before continuing!'
        exit()

    return    

def addStatisticalFluctuations(exp_histo):
    for bin_i in range(exp_histo.shape[0]):
        for bin_j in range(exp_histo.shape[1]):
            exp_histo[bin_i, bin_j] = np.random.poisson(exp_histo[bin_i, bin_j])
    return exp_histo
     
# compare whether two dictionaries are 'the same' #
def compare_dictionaries(dict1, dict2, verbose=False):

     # check stupid differences like not being of same type or same length or != None #
     if dict1 == None or dict2 == None:
         return False
     if type(dict1) is not dict or type(dict2) is not dict:
         return False

     shared_keys = set(dict2.keys()) & set(dict2.keys())
     if not ( len(shared_keys) == len(dict1.keys()) and len(shared_keys) == len(dict2.keys())):
         return False

    # ... if all these trivial things match, check all elements of dict #
     dicts_are_equal = True
     for key in dict1.keys():
         # if element itself is a dict again, do recursive call of compare_dictionaries() #
         if type(dict1[key]) is dict:
             dicts_are_equal = dicts_are_equal and compare_dictionaries(dict1[key],dict2[key])
         # ... otherwise try to somehow check whether these are the same for all possible types #
         else:
             tPrint = True
             # try a comparison - no matter what type the element is #
             try:
                 if type(dict1[key]) in (np.ndarray, np.array, list):
                     keys_fit          = np.all( [ np.all(dict1[key][i] == dict2[key][i]) for i in range(len(dict1[key])) ] )
                 else: keys_fit        = np.all(dict1[key] == dict2[key])
             # ... hard to catch all possible element types, so might not work in weird cases... #
             except:
                 keys_fit  = False
                 if verbose:
                     # the comparison should work for all 'standard' types as shown above, #
                     # if no standard type chosen, a comparison is difficult ...           #
                     print "Keys in dictionaries can not be compared: "
                     if verbose:
                         try:
                             print "\t-> Key: ", key
                             print "\t-> Type 1: ",type(dict1[key])
                             print "\t-> Type 2: ",type(dict2[key])
                             print "\t-> Value 1: ",dict1[key]
                             print "\t-> Value 2: ",dict2[key]
                         except:
                             print "\n"
                             print dict2.keys()
                             print dict1.keys()
                             tPrint = False
             dicts_are_equal = dicts_are_equal and keys_fit
             # if comparison failed for some reason, print reason (i.e. the key) why it failed, if verbose #
             if (verbose and not keys_fit) and tPrint: 
                 print "check failed: "+str(key)+": "+str(dict1[key])+ "(wanted) vs "+str(dict2[key]), "(found) -> ", dicts_are_equal

     return dicts_are_equal



# removes all keys from a dict-array-combination at certain depth (necessary to twist afterwards) #
def keepKeys(tDict, keys, depth):
    if depth == 0:
        if type(tDict) == dict:
            curKeys = tDict.keys()
            for tKey in curKeys:
                if not tKey in keys: del tDict[tKey]
        else:
            error( "Can not remove keys from non-dict type variable: "+str(type(tDict)) )
    else:
        for tKey in tDict:
            if type(tDict) == dict: keepKeys(tDict[tKey], keys, depth=depth-1)
            else:                   keepKeys(tKey, keys, depth=depth-1)


def compare_dictionaries_subset(dict1, dict2, verbose = False):
    if dict1 == None or dict2 == None:
        return False

    if type(dict1) is not dict or type(dict2) is not dict:
        return False

    shared_keys = set(dict2.keys()) & set(dict2.keys())
    
    dicts_are_equal = True
    for key in shared_keys:
        if type(dict1[key]) is dict:
            dicts_are_equal = dicts_are_equal and compare_dictionaries_subset(dict1[key],dict2[key])
        else:
            if verbose: print "check "+str(key)+": "+str(dict1[key])+ " vs "+str(dict2[key])
            if type(dict1[key]) == list:
                if len(dict1[key]) != len(dict2[key]):
                    return False
                for iList in range(len(dict1[key])):
                    dicts_are_equal = dicts_are_equal * (dict1[key][iList] == dict2[key][iList])
                    if (type(dicts_are_equal) == np.array or type(dicts_are_equal) == list or 
                        type(dicts_are_equal) == np.ndarray):
                        dicts_are_equal = np.prod(dicts_are_equal)
            else:
                dicts_are_equal = dicts_are_equal and (dict1[key] == dict2[key])
    return dicts_are_equal            
     

def getMatchingSet(dbFile, loaderDict, mcSettings, matchKeys, verbose):
    # Loop over the infoKey
    infile = open(dbFile)
    fileData   = pickle.load(infile)
    infile.close()

    floaderDict = deepcopy(loaderDict)
    fmcSettings = deepcopy(mcSettings)

    # Remove the information that will not be checked
    for oneKey in floaderDict.keys():
        if oneKey not in matchKeys:
            floaderDict.pop(oneKey)
    for oneKey in fmcSettings.keys():
        if oneKey not in matchKeys:
            fmcSettings.pop(oneKey)

    if verbose:
        print 'Testing for loader: ', floaderDict.keys()
        print 'Testing for mc settings:', fmcSettings.keys()

    # Loop over the settings stored
    for i in range(len(fileData['loader_settings'])):
        if (compare_dictionaries_subset(fileData['loader_settings'][i], floaderDict) and 
            compare_dictionaries_subset(fileData['mc_settings'][i], fmcSettings)):
            print 'Reading from file'
            return fileData['data'][i]

    # If I made it here, it didn't find anything suitable
    return False

    
def storeSet(dbFile, loaderDict, mcSettings, storeData):
    if os.path.exists(dbFile):
        infile = open(dbFile)
        fileData = pickle.load(infile)
        infile.close()
    else:
        fileData = {'loader_settings':[],
                    'mc_settings': [],
                    'data':     []}

    fileData['loader_settings'].append(loaderDict)
    fileData['mc_settings'].append(mcSettings)
    fileData['data'].append(storeData)
    
    pickle.dump(fileData, open(dbFile, 'w'))

    return 
