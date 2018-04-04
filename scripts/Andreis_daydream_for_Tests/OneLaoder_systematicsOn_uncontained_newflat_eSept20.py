

import os, sys
modules_dir = '/home/trwood/JP_fraction_original_jan26/modules'
#modules_dir = '/home/trwood/Andreis_daydream/modules'
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


#data_dir='/home/trwood/MSU_sample/MSU_sample_sept2016/data/'

#true_axis = np.array([[0., 12.11527659],[ 12.11527659,17.7827941],[17.7827941,26.10157216],[26.10157216,38.3118685],[38.3118685,56.23413252],[56.23413252,  np.inf]])


loader_nobkrd_a = dataLoader.dataLoader(observables =
                                ['reco_energy', 'reco_zenith', 'delta_llh'],
                                bin_edges   =
                                #no events below 10GeV so start binning at 1
                                #smarter way cut out the frist 3 bins.  
                                #retain the same binning but get rid of below 10GeV
                                #[10**np.linspace(0.75,1.75,9),
                                [10**np.linspace(0.5,2.25,12),
                                 np.arccos(np.linspace(-1.,1.,9))[::-1],
                                 np.array([-3, 2, np.inf])],	
				 #user = 'Chi2msu_no_background_noData_baselineONLY_flat',
                                 #user = 'Chi2msu_no_background_noData_systematicsON',
			#	user = 'TryOneLoadierSystemacis_flat_uncontained',
                                user = 'TryOneLoaderSystemacis_flat_uncontained',
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
                             
   #weight_keys = ['tweight_e', 'tweight_mu_k', 'tweight_mu_p'],
				#weight_keys = ['weight_e', 'weight_mu'], 
                                weight_keys = ['tweight_newflat_e', 'tweight_newflat_mu_k', 'tweight_newflat_mu_p'],
                                #extra_cuts = {'energy':true_axis[0],'bdt_score': (0.2,np.inf)},
		#		extra_cuts = {'bdt_score': (0.2,np.inf), 'mn_stop_contained': (0.1, 2.0) },
                                detsys_redo = False,
                                verbose = True,
                                table_nbins = 150)
