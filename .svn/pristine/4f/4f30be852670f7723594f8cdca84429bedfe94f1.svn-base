#!/usr/bin/env python


import os, sys
modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/modules'
sys.path.append(modules_dir)
modules_dir2 = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/user_settings'
sys.path.append(modules_dir2)

from numpy import *
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

    fit_results = ['llh', 'dm31','theta23','theta13','mix_angle','norm_nu','norm_e', 'norm_atmmu',
                   'atmmu_f','noise_f', 'norm', 'nu_nubar',
                   'norm_tau','norm_noise','nubar_ratio','uphor_ratio','gamma','axm_qe','axm_res',
                   'pid_bias','domeff','hole_ice','had_escale','norm_nc', 'chi_squared']

    aux_params = ['chi_squared','llh','norm_nu','norm_atmmu','norm_noise']

    # Initialize scan results
    scan_results = {}
    for fit_key in fit_results:
        scan_results[fit_key] = zeros([len(scan_info['theta23_list']), len(scan_info['dm31_list'])])

    missing_files = []
    nan_files     = []
    index = 0
    error = 0

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
                if one_fit['llh'] == None or isnan(one_fit['llh']):
                    print one_fit['llh']
                    nan_files.append(index+1)
                else:
                    for fit_key in fit_results:
                        scan_results[fit_key][theta23_i, dm31_i] = one_fit[fit_key]
            else:
                missing_files.append(index+1)

            index += 1

    print '\nLLH scan: Missing files -',len(missing_files), missing_files
    print 'LLH scan: NAN files -', len(nan_files), nan_files

    missing_files = missing_files+ nan_files

    if len(missing_files) > 0:
        print 'LLH scan: Trying to reprocess missing files with mode ', options.MODE, '. Continue?'
        raw_input()

        submit_string = ' '.join(['./run_scan.py -n', options.NAME, '-f 1 ', '-m',options.MODE,'-j ']) + \
            ','.join(["%i" % x for x in missing_files])
        print submit_string
        os.system(submit_string)
        exit()

    # Get the LLH miniumum here

    min_llh = scan_results['llh'].min()
    min_position = unravel_index(scan_results['llh'].argmin(), scan_results['llh'].shape)
    print '\n***Parameters at best fit***'
    scan_results['bestfit_scan'] = {} 
    for fit_key in fit_results:
        scan_results['bestfit_scan'][fit_key] = scan_results[fit_key][min_position]
        print fit_key, '\t', scan_results['bestfit_scan'][fit_key]
        # Also re-setting the fit_settings for the refined scan. Not passing the derived parameters.
        if fit_key in aux_params:
            continue
	try:
            scan_info['fit_settings'][fit_key][0] = scan_results['bestfit_scan'][fit_key]
	except:
    	    print 'shooot' 

    # # All results saved
    scan_results['dm31_list']    = scan_info['dm31_list']
    scan_results['theta23_list'] = scan_info['theta23_list']


    if options.REFIT:
        print '\nLLH scan: Getting the global fit again using the best fit from the scan as seed'
        import oscFit, dataLoader
        data_loaders = []
        for loader_settings in scan_info['loader_settings']:
            data_loaders.append(dataLoader.dataLoader(**loader_settings))
        fitter = oscFit.fitOscParams()
        overall_results = fitter(data_histograms = scan_info['data_histograms'],
                                 data_loaders    = data_loaders,
                                 fit_settings    = scan_info['fit_settings'],
                                 store_fit_details = True)

        scan_results['bestfit'] = overall_results
    else:
        print '\nLLH scan: Taking the global fit from the scan'
        scan_results['bestfit'] = scan_results['bestfit_scan']

    # Include everything in the scan results, including the settings.
    scan_results.update(scan_info)

    f = open(os.path.join(indir,'ScanResults.pckl'), 'w')
    pickle.dump(scan_results, f)
    f.close()


if __name__ == '__main__':
    from optparse import OptionParser


    usage = "usage: %prog [options] inputfile"
    parser = OptionParser(usage)
    parser.add_option("-n", "--name", default=None,
                      dest="NAME", help = "Name of the scan created")
    parser.add_option("-f", "--figures", action="store_true", dest='FIGURES',
                      help='Collect the results and produce figures')
    parser.add_option("-r", "--refit", action="store_true", dest="REFIT",
                      help="Re-do a global fit to refine the best fit point")
    parser.add_option("-m", "--mode", default="local",
                      dest="MODE", help = "Select where to run missing scan points")

    (options,args) = parser.parse_args()

    collectResults(options)
