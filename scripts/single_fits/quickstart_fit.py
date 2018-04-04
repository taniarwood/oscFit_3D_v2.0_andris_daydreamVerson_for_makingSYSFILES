#!/usr/bin/env python


import os, sys

# Adding the path of the modules
modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/modules'
sys.path.append(modules_dir)

# Importing modules
import pickle
import dataLoader, oscFit
import numpy as np


loader =  dataLoader.dataLoader(observables = 
                                ['reco_energy', 'reco_zenith', 'tracklength'],
                                bin_edges   =    
                                [10**np.linspace(0.8,1.75,9),
                                 np.arccos(np.linspace(-1,0.,9))[::-1],
                                 np.array([-np.inf, np.inf])],
                                user = 'micheal_use_pretania', 
                                LEaxis = [],#np.linspace(1, 3.2, 21),      
                                weight_keys = ['tweight_e', 'tweight_mu_k', 'tweight_mu_p'],
                                detsys_redo = True,
                                verbose = False,
				table_nbins = 150)

data_settings = {  'dm31':                     0.0025, 
                   'theta23':                  0.78,
                   'theta13':                  0.155,
                   'norm_nu':                  3.6, # In years!
                   'norm_e':                   1.,
                   'domeff':                   1., 
                   'nu_pi_scale':              0.6,
                   'hole_ice':                 0.02, 
                   'atmmu_template':           'data',
                   'simulation':               'baseline', 
		   'oscTables':		       'True',
                   'oscMode':                  'Vacuum',  
                   'ma_variations':            True,   
                   'add_detector_systematics': True}

data_histogram = loader.loadMCasData(data_settings,
                                     statistical_fluctuations=False)

print 'This is your data'
print data_histogram

fitter = oscFit.fitOscParams()
data_loaders    = [loader]          
data_histograms = [data_histogram]  

# The settings of the fit are given as [ini_value, fix?]
fit_settings    = {
    'simulation':     'baseline',      
    'dm31':           [0.0025, False, 'NH'], 
    'theta23':        [0.78, False],         
    'theta13':        [0.155, True],        
    'oscMode':        'TwoNeutrino',         
    'oscTables':      True,         #used to be False.  need to think about the tables.  need to understand resolution to set n bins
    'norm':           [1., True],
    'norm_e':         [1., False],
    'norm_tau':       [1., True],            
    'nu_nubar':       [1., True],
    'nubar_ratio':     [0., True],
    'uphor_ratio':    [0., True],
    'nu_pi_scale':    [1., False],

    'gamma':          [0., False],  
    'axm_qe':         [0., True],
    'axm_res':        [0., True],
    'pid_bias':       [0., True],
    'hole_ice':       [0.02, False], 

    'norm_nc':        [1., True],
    'domeff':         [1., False],
    'had_escale':     [1., True],
    'atmmu_f':        [0.05, False, 'mc'],
    'noise_f':        [0.0, True],
    'detector_syst':  True,
    'include_priors': True, 
    'printMode':      -1}                  

fit_priors = {'hole_ice':[0.02, 0.01],
              'gamma':[0.05, 0.1],
              'norm_e':[1., 0.2]}

# Obtain the result
result = fitter(data_histograms=data_histograms,
                data_loaders=data_loaders,
                fit_settings=fit_settings,
                fit_priors = fit_priors,
                ncalls = 1000,
                store_fit_details = True)
