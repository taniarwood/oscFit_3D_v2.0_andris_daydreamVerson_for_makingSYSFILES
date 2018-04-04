#!/usr/bin/env python

import os, sys
modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/modules'
sys.path.append(modules_dir)
modules_dir2 = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/user_settings'
sys.path.append(modules_dir2)

import numpy as np
import pickle
import oscFit_settings as user


############
## DATA AND FITTER SETTINGS
############
id_string = 'TestScan_MC_2PID_LE_Prob3'
oscMode   = 'Prob3'

outdir      = os.path.join(user.scans_dir, id_string)
outdir_bins = os.path.join(outdir, 'fits')
if os.path.isdir(outdir):
    print '\nLLH scan: Directory already exists. Do you want to remove the directory', id_string,'?'
    raw_input()
if not os.path.isdir(outdir):
    os.mkdir(outdir)
if not os.path.isdir(outdir_bins):
    os.mkdir(outdir_bins)


import dataLoader


loader =  dataLoader.dataLoader(bin_edges   = [10**np.linspace(0.8,1.75, 9),
                                               np.arccos(np.linspace(-1,0.,9))[::-1],
                                               np.array([-np.inf, 0.7, np.inf])],
                                observables = ['reco_energy', 'reco_zenith', 'pid'],
                                user = 'jpall_newdetector', 
                                LEaxis = [],
                                weight_keys = ['weight_e', 'weight_mu'], 
                                verbose = False,
                                detailed_detsys = False)


data_settings = {'dm31':                     0.0026, 
                 'theta23':                  0.8,
                 'theta13':                  0.155,
                 'mix_angle':                0.97, 
                 'norm_nu':                  2.4, # In years
                 'norm_e':                   1.,
                 'norm_tau':                 1.,
                 'norm_noise':               0.,
                 'noise_f':                  False, #Set either a normalization or fraction
                 'norm_nc':                  1.,
                 'norm_atmmu':               False, 
                 'atmmu_f':                  0.05, #Set either a normalization or fraction
                 'nu_nubar':                 1., #
                 'gamma':                    0.05, 
                 'axm_qe':                   0.,
                 'axm_res':                  0.,
                 'pid_bias':                 0.0,
                 'domeff':                   1., 
                 'hole_ice':                 0.02, 
                 'had_escale':               1.,
                 'atmmu_template':           'data',
                 'simulation':               'baseline', 
                 'oscMode':                  oscMode,
                 'oscTables':                False,
                 'ma_variations':            True,
                 'add_detector_systematics': True}


fit_settings  = {
    'simulation':     'baseline',
    'dm31':           [0.0027, False, 'NH'],
    'theta23':        [0.83, False], 
    'theta13':        [0.155, False],
    'mix_angle':      [1.0, False, 1.5],
    'oscMode':        oscMode,
    'oscTables':      False,
    'norm':           [1., False],
    'norm_e':         [1., False],
    'norm_tau':       [1., True], 
    'gamma':          [0.04, False],  
    'axm_qe':         [0., False],
    'axm_res':        [0., False],
    'pid_bias':       [0., True],
    'norm_nc':        [1., False],
    'hole_ice':       [0.02, False], 
    'domeff':         [1., False],
    'had_escale':     [1., True],
    'atmmu_f':        [0.07, False, 'data'],
    'noise_f':        [0.0, True],
    'fix_norm_region':[1.7],
    'detector_syst':  True,
    'include_priors': True, 
    'printMode':      True}



############
## PARAMETER SCAN SETTINGS
############

# If you are running in TwoNeutrino mode interpret these values as *mixing angle* (between 0 and 1)
sintheta23_2_i = 0.25
sintheta23_2_f = 0.75
sintheta23_2_points = 21

dm31_i = 1.8E-3
dm31_f = 3.6E-3
dm31_points = 21


theta23_list = np.arcsin(np.sqrt(np.linspace(sintheta23_2_i, sintheta23_2_f, sintheta23_2_points)))
dm31_list  = np.linspace(dm31_i, dm31_f, dm31_points)

dm31_map = np.repeat(dm31_list, len(theta23_list))
theta23_map = theta23_list.tolist()*len(dm31_list)

print 'LLH scan: Points in scan ', dm31_map.size


total_jobs = len(dm31_map)

scan_info = {'loader_settings':[loader.iniDict],
             'data_settings':data_settings,
             'fit_settings':fit_settings,
             'dm31_list':dm31_list, 'theta23_list':theta23_list,
             'dm31_map':dm31_map, 'theta23_map':theta23_map}

f = open(os.path.join(outdir,'OscFit_ScanSettings.pckl'), 'w')
pickle.dump(scan_info,f )
f.close()

print '\nLLH scan: Settings created!'
print '\tNow you can run the scan with run_mctest_scan.py -n ', id_string
print '\tSelect where to run the scan (-m) and if you want to run a test point first (-t)\n'
