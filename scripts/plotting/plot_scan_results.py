#!/usr/bin/env python


import os, sys
modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/modules'
sys.path.append(modules_dir)
modules_dir2 = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/user_settings'
sys.path.append(modules_dir2)

import oscFit_settings as user

import matplotlib
matplotlib.use('Agg')

import plotTools
import pickle


# You can pass a list of scans
# Figures are stored in the directory of the first scan listed
scan_names = sys.argv[1:]

scan_results = []

for one_scan in scan_names:
    scan_results.append(pickle.load(open(os.path.join(user.scans_dir, one_scan, 'ScanResults.pckl'))))

# You can re-define your labels and colors here
iclabels = scan_names
iccolor   = ['#283A90','b','g','b','k','m']

# Just the results given in scan names
plotTools.plotScanResults(scan_results, 
                          outdir = os.path.join(user.scans_dir, scan_names[0]), 
                          iclabels = iclabels, 
                          compare_results =  False, 
                          iccolor = iccolor)

# Compare against established experiments
plotTools.plotScanResults(scan_results, 
                          outdir = os.path.join(user.scans_dir, scan_names[0]), 
                          iclabels = iclabels, 
                          compare_results =  True, 
                          iccolor = iccolor)

print '\nData/MC comparisons for all', len(scan_results), 'scan(s) will be produced now.' 
print 'These comparisons take a couple of minutes to produce. Proceed?'
raw_input() 
# Produce the best fit point

scan_settings = pickle.load(open(os.path.join(user.scans_dir, one_scan, 'OscFit_ScanSettings.pckl')))

scan_results_forplots = {'result': scan_results[0]['bestfit'],
                         'loader_settings':scan_settings['loader_settings'],
                         'data_histograms':scan_settings['data_histograms'],
                         'data_settings':scan_settings['data_settings']}
plotTools.plotBestFit3D(scan_results_forplots, 
                        outdir = os.path.join(user.scans_dir, scan_names[0]))
