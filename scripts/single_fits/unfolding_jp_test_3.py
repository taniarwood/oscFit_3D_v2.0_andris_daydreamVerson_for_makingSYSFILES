

import os, sys
modules_dir = '/home/trwood/JP_fraction_original_jan26/modules'
sys.path.append(modules_dir)


# Importing modules
import argparse
import math
import pickle
import dataLoader, oscFit
import numpy as np
import oscfit_default_values as defs
from copy import deepcopy
import jp_mpl as jplot


data_dir='/home/trwood/MSU_sample/MSU_sample_sept2016/data/'

#true_axis = np.array([[5., 12.11527659],[ 12.11527659,17.7827941],[17.7827941,26.10157216],[26.10157216,38.3118685],[38.3118685,56.23413252],[56.23413252,  177.827941]])
#true_axis = np.array([[0., 12.11527659],[ 12.11527659,17.7827941],[17.7827941,26.10157216],[26.10157216,38.3118685],[38.3118685,56.23413252],[56.23413252,  np.inf]])
true_axis = np.array([[0, 8.],[8., 15.], [15., 25.],[25., 40.], [40., 80.], [80., 1000.]])

loader_nobkrd_a = dataLoader.dataLoader(observables =
                                ['reco_energy', 'reco_zenith', 'delta_llh'],
                                bin_edges   =
                                #no events below 10GeV so start binning at 1
                                #smarter way cut out the frist 3 bins.  
                                #retain the same binning but get rid of below 10GeV
                                  #was 10**np.linspace([0.75,1.25,2])
                                [10**np.linspace(0.75,2.25,11),
                                 np.arccos(np.linspace(-1.,1.,9))[::-1],
                                 np.array([-3, 2, np.inf])],	
				 user = 'Chi2msu_no_background_noData_baselineONLY_flat',
                                 #user = 'Chi2msu_no_background_noData_systematicsON',
                                LEaxis = [],#np.linspace(1, 3.2, 21),      
                                #tells the fit if you want to use the polynomial or not
                                #do linear instead (see one at end of brackets
                #ask about nuspecs and muspecs.. not in msu original user i have as reference
                #but i want this to be as similar to my larson evaluation as possible
                                detsys_nuspecs = {},
                                detsys_muspecs = {},
                                weight_keys = ['tweight_newflat_e', 'tweight_newflat_mu_k', 'tweight_newflat_mu_p'],
                                extra_cuts = {'energy':true_axis[0]},
                                detsys_redo = False,
                                verbose = True,
                                table_nbins = 150)



#make a fully new copy of loader_dict, ie don't just copy the location of the pointer
from copy import deepcopy
loader_dict = deepcopy(loader_nobkrd_a.iniDict)


loader_dict['extra_cuts'] =  {'energy':true_axis[1]}
loader_nobkrd_b= dataLoader.dataLoader(**loader_dict)


loader_dict['extra_cuts'] =  {'energy':true_axis[2]}
loader_nobkrd_c= dataLoader.dataLoader(**loader_dict)


loader_dict['extra_cuts'] =  {'energy':true_axis[3]}
loader_nobkrd_d= dataLoader.dataLoader(**loader_dict)


loader_dict['extra_cuts'] =  {'energy':true_axis[4]}
loader_nobkrd_e= dataLoader.dataLoader(**loader_dict)


loader_dict['extra_cuts'] =  {'energy':true_axis[5]}
loader_nobkrd_f= dataLoader.dataLoader(**loader_dict)




loader_dict['extra_cuts'] =  {}
#loader_dict['weight_keys'] = ['tweight_e', 'tweight_mu_k', 'tweight_mu_p']
loader_dict['user'] = 'Chi2noSystematicSets_new_finalized_msu_20122014_holeice2_baselineONLY_flat'
loader_bkrd_nu= dataLoader.dataLoader(**loader_dict)




loader_dict['extra_cuts'] =  {}
loader_dict['user'] = 'Chi2msu_BKGRND_only_flat'
loader_pureBkrd = dataLoader.dataLoader(**loader_dict)



#need to check that my fit settings that I am fixing match the data settings
fit_settings = deepcopy(defs.default_fit_settings)
fit_settings_fix    = {
    'simulation':     'baseline',
    'dm31':           [0.0025, False, 'NH'],
    'theta23':        [0.74, False],
    'theta13':        [0.155, False],
    'oscMode':        'TwoNeutrino',
    'oscTables':      False,         #used to be False.  need to think about the tables.  need to understand resolution to set n bins
    'norm':           [1., True],    #try once with this true and once with this false, same settings

    # These are good starting values for the flat spectrum
    #'nu_frac1':       [0.317, False],
    #'nu_frac2':       [0.168, False],
    #'nu_frac3':       [0.162, False],
    #'nu_frac4':         [0.138, False],
    #'nu_frac5':        [0.098, False],

    # These are good starting values for the realistic spectrum
    'nu_frac1':       [0.1614, False],
    'nu_frac2':       [0.26171, False],
    'nu_frac3':       [0.22058, False],
    'nu_frac4':         [0.16759, False],
    'nu_frac5':        [0.13429, False],



    'norm_e':         [1., False],
    'norm_tau':       [1., True],
    'nu_nubar':       [1., True],
    'nubar_ratio':    [0., True],
    'uphor_ratio':    [0., True],
    'nu_pi_scale':    [1., False],
    'hi_fwd':         [0.0, True],
    'gamma':          [0., True],
    'axm_qe':         [0., True],
    'axm_res':        [0., False],
    'pid_bias':       [0., True],
    'hole_ice':       [0.02, True],           # Fix these to true for this test???
    
    'mix_angle':      [1.0,False, 1.],

    'norm_nc':        [1., False],
    'domeff':         [1., True],
    'had_escale':     [1., True],
    'atmmu_f':        [0.0, True, 'data'],
    'noise_f':        [0.0, True],
    'detector_syst':  True,
    'include_priors': True,
    'printMode':      -1}  # (-1) this is printing each step


data_histo_bkrd_PlusNu = pickle.load(open('jp_data_real_spectrum.pckl'))
#data_histo_bkrd_PlusNu = pickle.load(open('jp_data.pckl'))


data_settings = deepcopy(defs.default_data_settings)
data_settings.update({'dm31':                  0.0025,
                   'theta23':                  0.74,
                   'theta13':                  0.155,
                   'mix_angle':                1.,
                   'norm_tau':                 1.,
                   'axm_qe':                   0.,
                   'axm_res':                  0.,
                   'norm_nu':                  2.4, # In years!
                   'norm_e':                   1.,
                   'domeff':                   1.,
                   'nu_pi_scale':              1.,

                   'hole_ice':                 0.02,
                   'atmmu_template':           'data',
                   'simulation':               'baseline',
                   'oscTables':                False,
                   'gamma':                    0.,
                   'ma_variations':            False,
                   'add_detector_systematics': False,
                   'norm_atmmu':               0.,
                   'pid_bias':                 0.,
                   'domeff':                   1.,
                   'hole_ice':                 0.02,
                   'hi_fwd':                   0., # MSU parameterization of the forward impact angle
                   'had_escale':               1.,
                   'oscMode':                  'TwoNeutrino', 
                   'ma_variations':            False,  
                   'norm_noise':                0.0,
                    'atmmu_f':                  0.0   })



d1 = loader_nobkrd_a.loadMCasData(data_settings, statistical_fluctuations=False)
d2 = loader_nobkrd_b.loadMCasData(data_settings, statistical_fluctuations=False)
d3 = loader_nobkrd_c.loadMCasData(data_settings, statistical_fluctuations=False)
d4 = loader_nobkrd_d.loadMCasData(data_settings, statistical_fluctuations=False)
d5 = loader_nobkrd_e.loadMCasData(data_settings, statistical_fluctuations=False)
d6 = loader_nobkrd_f.loadMCasData(data_settings, statistical_fluctuations=False)
bkg = loader_pureBkrd.loadMCasData(data_settings, statistical_fluctuations=False)

print '\nCHECKING THE NUMBER OF EVENTS IN DATA'
print np.sum(data_histo_bkrd_PlusNu)
print 'This is the sum of all the individual histograms'


all_data = d1+d2+d3+d4+d5+d6
print np.sum(all_data)

print '\n THESE ARE THE INDIVIDUAL NUMBERS'
print d1
print '\n'
print d2 
print '\n'
print d3
print '\n'
print d4
print '\n'
print d5
print '\n'
print d6






import iminuit
fitter = oscFit.fitOscParams()


fit_priors = {'hole_ice':[0.02, 0.01],
  #            'gamma':[0.05, 0.1],
              'norm_e':[1., 0.2]}


oscFitResults_fixed =fitter(data_histograms=[data_histo_bkrd_PlusNu],
                            data_loaders=[loader_nobkrd_a,
                                          loader_nobkrd_b,
                                          loader_nobkrd_c,
                                          loader_nobkrd_d,
                                          loader_nobkrd_e,
                                          loader_nobkrd_f,
                                          loader_pureBkrd],
                            fit_settings=fit_settings_fix,
                            ncalls               = 1000, # For migrad. Simplex will use twice the number.
                            store_fit_details    = True)
#  return_minuit_object = False) 

pickle.dump(oscFitResults_fixed,open('/home/trwood/oscfit_output/msu_repickled_tania/JP_test_3postMadison.pckl', 'w'))

#pickle.dump(oscFitResults_fixed,open('/home/trwood/oscfit_output/msu_repickled_tania/JP_test_3'  + str(norm_1_use) +'postMadison.pckl', 'w'))

#np.savetxt('/home/trwood/oscfit_output/msu_repickled_tania/JP_test_3.' + str(norm_1_use) +'.postMadison.txt',
           np.c_[float(oscFitResults_fixed['nu_pi_scale']), oscFitResults_fixed['llh']], fmt='%1.4f', delimiter=' ',
           newline=os.linesep)
