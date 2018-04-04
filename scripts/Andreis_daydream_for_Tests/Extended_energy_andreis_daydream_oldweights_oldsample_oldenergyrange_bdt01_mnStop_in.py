

import os, sys
#modules_dir = '/home/trwood/JP_fraction_original_jan26/modules'
modules_dir = '/home/trwood/Andreis_daydream/modules'
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

#true_axis = np.array([[0., 12.11527659],[ 12.11527659,17.7827941],[17.7827941,26.10157216],[26.10157216,38.3118685],[38.3118685,56.23413252],[56.23413252,  np.inf]])


loader_nobkrd_a = dataLoader.dataLoader(observables =
                                ['reco_energy', 'reco_zenith', 'delta_llh'],
                                bin_edges   =
                                #no events below 10GeV so start binning at 1
                                #smarter way cut out the frist 3 bins.  
                                #retain the same binning but get rid of below 10GeV
                               # [10**np.linspace(0.75,1.75,9),
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
                                #weight_keys = ['tweight_e', 'tweight_mu_k', 'tweight_mu_p'],
				weight_keys = ['weight_e', 'weight_mu'], 
                                #weight_keys = ['tweight_newflat_e', 'tweight_newflat_mu_k', 'tweight_newflat_mu_p'],
                                #extra_cuts = {'energy':true_axis[0],'bdt_score': (0.2,np.inf)},
				extra_cuts = {'bdt_score': (0.1,np.inf)},
                                detsys_redo = False,
                                verbose = True,
                                table_nbins = 150)



#make a fully new copy of loader_dict, ie don't just copy the location of the pointer
from copy import deepcopy
loader_dict = deepcopy(loader_nobkrd_a.iniDict)

'''
loader_dict['extra_cuts'] =  {'energy':true_axis[1],'bdt_score': (0.2,np.inf)}
loader_nobkrd_b= dataLoader.dataLoader(**loader_dict)


loader_dict['extra_cuts'] =  {'energy':true_axis[2],'bdt_score': (0.2,np.inf)}
loader_nobkrd_c= dataLoader.dataLoader(**loader_dict)


loader_dict['extra_cuts'] =  {'energy':true_axis[3],'bdt_score': (0.2,np.inf)}
loader_nobkrd_d= dataLoader.dataLoader(**loader_dict)


loader_dict['extra_cuts'] =  {'energy':true_axis[4],'bdt_score': (0.2,np.inf)}
loader_nobkrd_e= dataLoader.dataLoader(**loader_dict)


loader_dict['extra_cuts'] =  {'energy':true_axis[5],'bdt_score': (0.2,np.inf)}
loader_nobkrd_f= dataLoader.dataLoader(**loader_dict)




#loader_dict['extra_cuts'] =  {}
##loader_dict['weight_keys'] = ['tweight_e', 'tweight_mu_k', 'tweight_mu_p']
#loader_dict['user'] = 'Chi2noSystematicSets_new_finalized_msu_20122014_holeice2_baselineONLY_flat'
#loader_bkrd_nu= dataLoader.dataLoader(**loader_dict)




loader_dict['extra_cuts'] =  {}
loader_dict['user'] = 'Chi2msu_BKGRND_only_flat'
#loader_pureBkrd = dataLoader.dataLoader(**loader_dict)

'''
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
  #                 'nu_pi_scale':              1.,

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




#need to check that my fit settings that I am fixing match the data settings
fit_settings = deepcopy(defs.default_fit_settings)
fit_settings_fix    = {
    'simulation':     'baseline',
    'dm31':           [0.0025, True, 'NH'],
    'theta23':        [0.74, False],
    'theta13':        [0.155, False],
    'oscMode':        'TwoNeutrino',
    'oscTables':      False,         #used to be False.  need to think about the tables.  need to understand resolution to set n bins
    'norm':           [1., True],    #try once with this true and once with this false, same settings
#    'nu_frac1':       [0.3, False],
#    'nu_frac2':       [0.2, False],
#    'nu_frac3':       [0.2, False],
#    'nu_frac4':         [0.15, False],
#    'nu_frac5':        [0.1, False],
#  #  'nu_frac6':        [0.1, False],
#  #  'nu_frac7':        [0.1, False],
#  #  'nu_frac8':        [0.1, False],
#  #  'nu_frac9':        [0.1, False],
    'norm_e':         [1., True],
    'norm_tau':       [1., True],
    'nu_nubar':       [1., True],
    'nubar_ratio':    [0., True],
    'uphor_ratio':    [0., True],
 #   'nu_pi_scale':    [1., True],
    'hi_fwd':         [0.0, True],
    'gamma':          [0., True],
    'axm_qe':         [0., True],
    'axm_res':        [0., True],
    'pid_bias':       [0., True],
    'hole_ice':       [0.02, True],           # Fix these to true for this test???
    
    'mix_angle':      [1.0,True, 1.],

    'norm_nc':        [1., True],
    'domeff':         [1., True],
    'had_escale':     [1., True],
    'atmmu_f':        [0.0, True, 'data'],
    'noise_f':        [0.0, True],
    'detector_syst':  True,
    'include_priors': True,
    'printMode':      -1}  # (-1) this is printing each step


#data_histo_bkrd_PlusNu =loader_bkrd_nu.loadMCasData(data_settings,
#                                         statistical_fluctuations=False)

d1 = loader_nobkrd_a.loadMCasData(data_settings, statistical_fluctuations=False)
#d2 = loader_nobkrd_b.loadMCasData(data_settings, statistical_fluctuations=False)
#d3 = loader_nobkrd_c.loadMCasData(data_settings, statistical_fluctuations=False)
#d4 = loader_nobkrd_d.loadMCasData(data_settings, statistical_fluctuations=False)
#d5 = loader_nobkrd_e.loadMCasData(data_settings, statistical_fluctuations=False)
#d6 = loader_nobkrd_f.loadMCasData(data_settings, statistical_fluctuations=False)
#bkg = loader_pureBkrd.loadMCasData(data_settings, statistical_fluctuations=False)

print '\nCHECKING THE NUMBER OF EVENTS IN DATA'
#print np.sum(data_histo_bkrd_PlusNu)
print 'This is the sum of all the individual histograms'

all_data = d1
#all_data = d1+d2+d3+d4+d5+d6
print np.sum(all_data)


import pickle
pickle.dump(all_data,open('daydream_Joweight_data_real_spectrum.pckl', 'w'))
#pickle.dump(all_data,open('jp_data.pckl', 'w'))


print 'Done'
sys.exit()



import iminuit
fitter = oscFit.fitOscParams()


fit_priors = {'hole_ice':[0.02, 0.01],
  #            'gamma':[0.05, 0.1],
              'norm_e':[1., 0.2]}


oscFitResults_fixed, Minuit2Results_fixed=fitter(data_histograms=[data_histo_bkrd_PlusNu],
                                                 data_loaders=[loader_nobkrd_a,loader_nobkrd_b,loader_nobkrd_c,loader_nobkrd_d,loader_nobkrd_e,loader_nobkrd_f,loader_pureBkrd],
                                                 fit_settings=fit_settings_fix,
                                                 ncalls               = 1000, # For migrad. Simplex will use twice the number.
                                                 store_fit_details    = True)
#  return_minuit_object = False) 

oscFitResults_fixed.update({'data_settings': data_settings})

pickle.dump(oscFitResults_fixed,open('/home/trwood/oscfit_output/msu_repickled_tania/NoSytematics_may1_8bins_JasperTryGlobal_DRAGONChi'  + str(norm_1_use) +'_norm1234FreeNorm1Fixed.pckl', 'w'))

np.savetxt('/home/trwood/oscfit_output/msu_repickled_tania/NoSytematics_may1_8bins_JasperTryGlobal_Scan.' + str(norm_1_use) +'.oscFit_chi2_Results_nu_pi_scale.txt',
           np.c_[float(oscFitResults_fixed['nu_pi_scale']), oscFitResults_fixed['llh']], fmt='%1.4f', delimiter=' ',
           newline=os.linesep)
