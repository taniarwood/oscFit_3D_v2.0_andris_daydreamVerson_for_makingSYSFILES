

import os, sys
#modules_dir = '/home/trwood/JP_fraction_original_jan26/modules'
modules_dir = '/gs/project/ngw-282-ac/trwood/jasper_home/JP_fraction_original_jan26/modules/'
######## ie 'oscFit_3D_v2.0_andris_daydreamVerson_for_makingSYSFILES'= JP_fraction_original_jan26/, which is just oscFit3D_v2.0, checked out Jan 26th 2017
#modules_dir = '/home/trwood/Andreis_daydream/modules'
sys.path.append(modules_dir)

################################################################
#  Script that makes the SYSFILE for the H3a DPMJET-III version of the fitting machinery 
#march 2018

# Importing modules
import argparse
import math
import pickle
import dataLoader_4weights as dataLoader
import oscFit
import numpy as np
import oscfit_default_values as defs
from copy import deepcopy
import jp_mpl as jplot
import argparse

parser = argparse.ArgumentParser(description='script for calculating neutrino flux forward fold bin nomrmalization profiles')
#parser.add_argument('normbinfix', help='value to try the scan on')
parser.add_argument('crpmodel', help='cosmic ray primary model')
parser.add_argument('hintmodel', help='hadronic interaction model')

args = parser.parse_args()
pmodelu = args.crpmodel
intmodelu = args.hintmodel   #_model_tests


name = pmodelu + '_' + intmodelu
name_e_k =  name + '_e_k'
name_e_p =  name + '_e_p'
name_mu_k = name + '_mu_k'
name_mu_p = name + '_mu_p'

#USER FOR OSFIT dataloader:
user_name = pmodelu+"_"+intmodelu+"_user_setting"

#change interaction model name to have no periods in it which 
print intmodelu

import re
intmodelu = re.sub('[.]', '_', intmodelu)

print intmodelu

sysfile_name = pmodelu+"_"+intmodelu+"_sysfile.pckl"
user_name_use =  pmodelu+"_"+intmodelu+"_user_setting"
print sysfile_name
#one loader to rule them all:
sysfile_name_use = sysfile_name 
loader_nobkrd_a = dataLoader.dataLoader(observables =
                                ['reco_energy', 'reco_zenith', 'delta_llh'],
                                bin_edges   =
                                #no events below 10GeV so start binning at 1
                                #smarter way cut out the frist 3 bins.  
                                #retain the same binning but get rid of below 10GeV
                                #[10**np.linspace(0.75,1.75,9),
                                [10**np.linspace(0.75,2.25,11),
                                 np.arccos(np.linspace(-1.,1.,9))[::-1],
                                 np.array([-3, 2, np.inf])],	
		                 user = user_name_use,
#                                user = 'TryOneLoaderSystemacis_flat_uncontained',
				#user = 'Chi2msu_no_background_noData_baselineONLY_flat_uncontained_numu_nue_both_albrecht_muons_too_DMP_h3a',
				#user = 'TryOneLoaderSystemacis_flat_uncontained_SYB23_GH',


				LEaxis = [],#np.linspace(1, 3.2, 21),      
                                #tells the fit if you want to use the polynomial or not
                                #do linear instead (see one at end of brackets
                #ask about nuspecs and muspecs.. not in msu original user i have as reference
                #but i want this to be as similar to my larson evaluation as possible
                            #    detsys_nuspecs = {},
                            #    detsys_muspecs = {},
				 detsys_nuspecs = {'domeff': [1., 1],
                                                     'hole_ice': [0.02, 1],
                                                     'hi_fwd': [0.0, 1]},     #[central value, polynomial degree]
                                detsys_muspecs = {},
                               # pmodel = pmodelu,
                               # crmodel = intmodelu,                           
   #weight_keys = ['tweight_e', 'tweight_mu_k', 'tweight_mu_p'],
				#weight_keys = ['weight_e', 'weight_mu'], 
                                #weight_keys = ['tweight_DMP_GH_e_jaspert', 'tweight_DMP_GH_mu_k_jaspert', 'tweight_DMP_GH_mu_p_jaspert'],
				#weight_keys = ['tweight_DMP_h3a_e_jaspert', 'tweight_DMP_h3a_mu_k_jaspert', 'tweight_DMP_h3a_mu_p_jaspert'],
				weight_keys = [name_e_k, name_e_p, name_mu_k, name_mu_p],
				#extra_cuts = {'bdt_score': (0.2,np.inf)},
		#		extra_cuts = {'bdt_score': (0.2,np.inf), 'mn_stop_contained': (0.1, 2.0) },
                                detsys_redo = False,
                                verbose = True,
                                table_nbins = -1,  #) #150)
				sysfile = sysfile_name_use)
