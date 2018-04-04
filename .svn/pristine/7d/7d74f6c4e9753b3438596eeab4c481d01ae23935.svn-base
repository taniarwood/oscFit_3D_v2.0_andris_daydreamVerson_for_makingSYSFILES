#!/usr/bin/env python


import os, sys
modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/modules'
sys.path.append(modules_dir)
modules_dir2 = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/user_settings'
sys.path.append(modules_dir2)

import matplotlib
matplotlib.use('Agg')

from numpy import *
import numpy as np
import jp_mpl as jplot
import matplotlib.pyplot as plt
import sys, os
import pickle
import oscFit_settings as user

def collectResults(options):
    indir = os.path.join(user.scans_dir, options.NAME)
    scan_info = pickle.load(open(os.path.join(indir,'OscFit_ScanSettings.pckl')))

    total_jobs = len(scan_info['theta23_map'])
    
    dm2_list     = scan_info['dm31_list']
    theta23_list = scan_info['theta23_list']



    fits_dir    = os.path.join(indir,'fits')
    infile_list = os.listdir(fits_dir)
    infile_list.sort()
    if len(infile_list) == 0:
        print '\n **** LLH scan: No files were found! Exiting ****'
        exit()


    fit_results = ['dm31','theta23','theta13','mix_angle','norm_e',
                   'atmmu_f','noise_f', 'norm_nu',
                   'norm_tau','numubar_e2','numubar_z2','gamma','axm_qe','axm_res',
                   'pid_bias','domeff','hole_ice','had_escale','norm_nc']

    skip_params = ['theta23', 'dm31','mix_angle']

    # Initialize scan results
    scan_results = {}
    for fit_key in fit_results:
        scan_results[fit_key] = []
    theta23_injected = []
    dm31_injected    = []

    missing_files = []
    nan_files     = []
    index = 0
    error = 0

    counter = 0
    # Lazy coding. Just look for the number
    for dm31_i in range(0, len(scan_info['dm31_list'])):
        if len(infile_list) == 0:
            break
        for theta23_i in range(0, len(scan_info['theta23_list'])):
            if len(infile_list) == 0:
                break
            required_task = 'Task_' + "%06i" % index + '.pckl'

            if required_task in infile_list:
                one_fit = pickle.load(open(os.path.join(fits_dir, required_task)))
                infile_list.pop(infile_list.index(required_task))
                if one_fit['llh'] == None or isnan(one_fit['llh']):# or not one_fit['successful_fit']:
                    #print one_fit['llh']
                    nan_files.append(index+1)
                else:
                    for fit_key in fit_results:
                        scan_results[fit_key].append( one_fit[fit_key])
                    dm31_injected.append(dm2_list[dm31_i])
                    theta23_injected.append(theta23_list[theta23_i])
            else:
                missing_files.append(index+1)

            index += 1

    theta23_injected = np.array(theta23_injected)
    dm31_injected = np.array(dm31_injected)
    print theta23_injected.shape
        
    print '\nLLH scan: Missing files -',len(missing_files)
    print 'LLH scan: NAN files -', len(nan_files)

    missing_files = missing_files+ nan_files

    if len(missing_files) > 0 and options.DOMISSING:
        print 'LLH scan: Trying to reprocess missing files with mode ', options.MODE, '. Continue?'
        raw_input()

        submit_string = ' '.join(['./run_mctest_scan.py -n', options.NAME,'-m',options.MODE,'-j ']) + \
            ','.join(["%i" % x for x in missing_files])
        print submit_string
        os.system(submit_string)
        exit()

    # Produce histograms and figures here

    # Theta23-Dm31 map
    fig = plt.figure()
    if scan_info['data_settings']['oscMode'] == 'TwoNeutrino':
        t23diff = np.array(theta23_injected) - np.array(scan_results['mix_angle'])
        plt.plot(theta23_injected, dm31_injected*10**3, 'xr')
        plt.plot(scan_results['mix_angle'], np.array(scan_results['dm31'])*10**3, '.b')
        plt.xlabel(r'$\sin^2(2\theta_{23})$', fontsize = 'large')
        t23key = 'Mixing angle'

    else:
        t23diff = np.sin(theta23_injected)**2 - np.sin(scan_results['theta23'])**2
        plt.plot(np.sin(theta23_injected)**2, dm31_injected*10**3, 'xr')
        plt.plot(np.sin(scan_results['theta23'])**2, np.array(scan_results['dm31'])*10**3, '.b')
        plt.xlabel(r'$\sin^2(\theta_{23})$', fontsize = 'large')
        t23key = 'sin(theta)**2'
    fig.savefig(os.path.join(indir,'FullScan.png'))

    fig = plt.figure()
    xaxis = np.linspace(np.min(t23diff), np.max(t23diff), 22)*1.1
    n, x = np.histogram(t23diff, xaxis)
    jplot.unfilledBar(x, n)
    plt.ylabel('Counts')
    plt.xlabel(t23key)
    fig.savefig(os.path.join(indir, 'Theta23.png'))

    fig = plt.figure()
    dm31diff = (dm31_injected - np.array(scan_results['dm31']))*10**3
    xaxis = np.linspace(np.min(dm31diff), np.max(dm31diff), 22)*1.1
    n, x = np.histogram(dm31diff, xaxis)
    jplot.unfilledBar(x, n)
    plt.ylabel('Counts')
    plt.xlabel('DM31 (10**3 eV)')
    fig.savefig(os.path.join(indir, 'DM31.png'))
        
               


    # Parameter-wise figures
    for one_key in fit_results:
        if one_key in skip_params:
            continue
        pdiff = np.array(scan_results[one_key]) - scan_info['data_settings'][one_key]
        if np.mean(pdiff) == 0:
            continue
        xaxis = np.linspace(np.min(pdiff), np.max(pdiff), 22)*1.1
        n, x = np.histogram(pdiff, xaxis)

        fig = plt.figure()
        jplot.unfilledBar(x, n)
        plt.xlabel(one_key)
        plt.ylabel('Counts')
        fig.savefig(os.path.join(indir, one_key + '.png'))




if __name__ == '__main__':
    from optparse import OptionParser


    usage = "usage: %prog [options] inputfile"
    parser = OptionParser(usage)
    parser.add_option("-n", "--name", default=None,
                      dest="NAME", help = "Name of the scan created")
    parser.add_option("-d", "--domissing", action="store_true", dest='DOMISSING',
                      help='Produce missing points before continuing')
    parser.add_option("-m", "--mode", default="local",
                      dest="MODE", help = "Select where to run missing scan points")

    (options,args) = parser.parse_args()

    collectResults(options)
