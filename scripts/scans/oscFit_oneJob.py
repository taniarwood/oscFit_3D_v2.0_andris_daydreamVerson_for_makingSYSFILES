import os, sys
modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/modules'
sys.path.append(modules_dir)
modules_dir2 = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/user_settings'
sys.path.append(modules_dir2)

import dataLoader, oscFit
import pickle
import oscFit_settings as user
import numpy as np

task_start      = int(sys.argv[1])
tasks_required  = int(sys.argv[2])
indir   = sys.argv[3]
indir = os.path.join(user.scans_dir, indir)

f = open(os.path.join(indir,'OscFit_ScanSettings.pckl'))
scan_settings = pickle.load(f)
f.close()

for task_nr in range(task_start*tasks_required, task_start*tasks_required + tasks_required):
    
    if task_nr >= len(scan_settings['dm31_map']):
        print 'Last job, finished!'
        exit()

    print '\nLLH scan: Doing task ', task_nr

    scan_settings['fit_settings']['dm31']      = [scan_settings['dm31_map'][task_nr], 
                                                  True, 
                                                  scan_settings['fit_settings']['dm31'][2]]
    
    if scan_settings['fit_settings']['oscMode'] == 'TwoNeutrino':
        scan_settings['fit_settings']['mix_angle'] = [scan_settings['theta23_map'][task_nr], True]
    else:
        scan_settings['fit_settings']['theta23']   = [scan_settings['theta23_map'][task_nr], True]

    # Need to obtain the data loaders
    data_loaders = []
    for loader_settings in scan_settings['loader_settings']:
        data_loaders.append(dataLoader.dataLoader(**loader_settings))

    if not scan_settings.has_key('fit_priors'):
        scan_settings['fit_priors'] = {}

    fitter = oscFit.fitOscParams()
    bestfit = fitter(data_histograms = scan_settings['data_histograms'],
                     data_loaders    = data_loaders,
                     fit_settings    = scan_settings['fit_settings'],
                     fit_priors      = scan_settings['fit_priors'],
                     store_fit_details = False) # Avoid storing correlation matrices for scans

    pickle.dump(bestfit, open(os.path.join(indir,'fits', 'Task_' + "%06i" % task_nr + '.pckl'), 'w'))
