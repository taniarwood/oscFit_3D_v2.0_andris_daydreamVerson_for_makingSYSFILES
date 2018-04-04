#!/usr/bin/env python


import os, sys

# Add the path where the modules of oscFit are located
modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/modules'
sys.path.append(modules_dir)

# Import some stuff necessary for declaring your things
import pickle
import dataLoader, oscFit
import numpy as np

# Store the results in a pickle file
store_result = False
if len(sys.argv) > 1:
    store_result= True
    result_name = sys.argv[1]
    outdir = '/afs/ifh.de/user/y/yanezjua/scratch/analysis-results/oscFit_single_results/'+result_name + '.pckl'
    print 'Storing the result as ', result_name
else:
    print 'Provide a name to store the result'

# Initialize instances of dataLoader. One per sample for the fit
# Note that oscFit3D already works in N-Dimensional histograms. 
# In principle one sample should be enough to include also a particle ID
# You can use multiple instances for fitting different datasets (e.g. IC79 + IC86)

loader =  dataLoader.dataLoader(observables =  # Define the name of the observables in the pickle file
                                ['reco_energy', 'reco_zenith', 'pid'],
                                bin_edges   =    #Define the edges that you want to use per observable
                                [10**np.linspace(0.8,1.75,9),
                                 np.arccos(np.linspace(-1,0.,9))[::-1],
                                 np.array([-np.inf, 0.7, np.inf])],
                                user = 'jpall_fast', # Name of the user
                                LEaxis = [],      # LE axis? To use LE, respect the order of observables given
                                weight_keys = ['weight_e', 'weight_mu'], # Name of the weight keys [nu_e, nu_mu]
                                extra_cuts = {},    # Data missing some cuts? Use 'variable':[llimit,hlimit]
                                break_energy = 1000.,# Where should GENIE stop and NuGen start? (GeV)
                                detsys_perflavor=False, # Calculate the detector variations per flavor
                                legacy_detsys = False, # Calcualte the detector variations as in the PRD paper
                                detailed_detsys=False, # Re-calculate the detector variations at each minimization step
                                verbose = False)    # Print lots of stuff?

# If you want to reproduce the same loader afterwards you can use dataLoader.dataLoader(**loader_loader.iniDict)
# The settings are saved to the pickle file. This is very effective for doing LLH scans in parallel.

# dataLoader takes care of putting all the MC files in shape for use.
# Now you need to produce your pseudo-data from MC, or get the actual data

# These settings are necessary for producing MC pseudo-data.
# The module dataLoader has default values already set.
# Look in there to have an idea what you can vary.
data_settings = {  'dm31':                     0.0025, 
              'theta23':                  0.78,
              'theta13':                  0.155,
              'mix_angle':                1., 
              'norm_nu':                  3., # In years!
              'norm_e':                   1.,
              'norm_tau':                 1.,
              'nu_nubar':                 1., #
              'nubar_ratio':              0.,
              'uphor_ratio':              0.,
              'gamma':                    0., 
              'axm_qe':                   0.,
              'axm_res':                  0.,
              'pid_bias':                 0.,
              'domeff':                   1., 
              'hole_ice':                 0.02, 
              'hi_fwd':                   0., # MSU parameterization of the forward impact angle
              'had_escale':               1.,
              'norm_nc':                  1.,
              'atmmu_template':           'data',
              'norm_atmmu':               0.01, 
              'simulation':               'baseline', 
              'oscMode':                  'Vacuum',  # Vacuum works by default.
              'ma_variations':            True,      # False if axial mass syst. are missing in pickle files
              'add_detector_systematics': False}

# This command will provide the MC data for an Asimov test. Data will look perfect.
# Setting statistical_fluctuations to True will add variations in the data. More realistic.
# Check that you are obtaining the number of events desired.
# Don't re-scale the data if you have set statistical fluctuations to True
data_histogram = loader.loadMCasData(data_settings,statistical_fluctuations=False)

print 'This is your data'
print data_histogram

# If you want to use the data instead, use this command
# data_settings = [] # Declare it empty if using real data
# data_histogram = loader.loadData(year = 'all')


# Now start an instance of the fitter oscFit

fitter = oscFit.fitOscParams()
data_loaders    = [loader]          # Here you could list many samples to be fit at once
data_histograms = [data_histogram]  # The data histogram for each of the samples

# The settings of the fit are given as [ini_value, fix?]
fit_settings    = {
    'simulation':     'baseline',      # This patches GENIE simulation with NuGen.
    'dm31':           [0.0025, False, 'NH'], # Need to specify the hierarchy of the fit
    'theta23':        [0.78, False],         # Not used in the "TwoNeutrino" mode
    'theta13':        [0.155, True],        # Not used in the "TwoNeutrino" mode
    'mix_angle':      [0.95, False, 1.4],     # Only used in the "TwoNeutrino" mode. Leave as is otherwise.
    'oscMode':        'TwoNeutrino',         # How do you want to fit the data?
    'oscTables':      True,
    'norm':           [1., True],
    'norm_e':         [1., False],
    'norm_tau':       [1., True],            # Set to false for nutau apperance studies
    'nu_nubar':       [1., True],
    'nubar_ratio':     [0., True],
    'uphor_ratio':    [0., True],

    'gamma':          [0., False],  
    'axm_qe':         [0., True],
    'axm_res':        [0., True],
    'pid_bias':       [0., True],
    'hole_ice':       [0.02, True], 
    'hi_fwd':         [0.0, True],

    'norm_nc':        [1., True],
    'domeff':         [1., True],
    'had_escale':     [1., True],
    'atmmu_f':        [0.05, False, 'data'],
    'noise_f':        [0.0, True],
    'fix_norm_region':[1.7],
    'detector_syst':  True,
    'include_priors': True, 
    'printMode':      True}                  # Do you want to see every step MINUIT takes?

# Define your own priors (optional)
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

# Save the result - All this info is necessary for plotting scripts!
# if store_result:
#     result_dict = {'result':result,
#                    'fit_settings':fit_settings,
#                    'loader_settings':[x.iniDict for x in data_loaders],
#                    'data_histograms':data_histograms,
#                    'data_settings':data_settings}
#     pickle.dump(result_dict, open(outdir, 'w'))
