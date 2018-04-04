import numpy as np
import sys, os
modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/modules'
sys.path.append(modules_dir)
print "Modules at: ", modules_dir
import pickle
theta24_inj = np.arcsin(np.sqrt(0.0))
theta34_inj = np.arcsin(np.sqrt(0.0))
mode = "Vacuum"

## Time: this script need ~100-110 seconds in Vacuum and ~500seconds in Prob3 mode


############
## FITTER SETTINGS
############
import dataLoader

energy_axis_edges = 10**np.linspace(0.8,1.75,9)
coszen_axis_edges = np.arccos(np.linspace(-1.,0.,9))[::-1]

tracks_loader =  dataLoader.dataLoader(
                observables = ['reco_energy',  'reco_zenith'],
                bin_edges   = [energy_axis_edges, coszen_axis_edges],
                user = 'PRD_new_lowdt_data', # Name of the user #PRD_extended_lowdt_data
                LEaxis = [],      # LE axis? To use LE, respect the order of observables given
                table_nbins  = 0,
                weight_keys = ['weight_e', 'weight_mu'], 
                extra_cuts = {
                              'reco_energy':[energy_axis_edges[0]-1e-5,
                                             energy_axis_edges[-1]+1e-5], 
                              'reco_zenith':[coszen_axis_edges[0]-1e-5,
                                             coszen_axis_edges[-1]+1e-5]
                              },    # Data missing some cuts? Use 'variable':[llimit,hlimit]
                break_energy = 1000.,
                verbose = True, 
                use_kde_bg       = False,
                use_kde_sys      = False,
                bootstrap_kde_bg = False,
                ) 



mc_settings = {'dm31':                  0.00258, 
            'theta23':                  0.82, 
            'mix_angle':                1.,
            'norm_e':                   1., 
            'norm_nu':                  2.6,
            'gamma':                    +0.01, 
            'domeff':                   1.0, 
            'hole_ice':                 0.02, 
            'norm_nc':                  1.,
            'atmmu_template':           'data',
            'norm_atmmu':               0.000,
            'simulation':               'baseline', 
            'oscMode':                  mode,
            'add_detector_systematics': True, 
            'theta24':                   0.0,
            'theta34':                   0.0,
            'dm41':                      1.0 }

tracks_data_histogram   = tracks_loader.loadMCasData(mc_settings,statistical_fluctuations=False)


#############
## MINIMIZATION
#############
fit_settings    = {            
            'simulation':     'baseline',
            'dm31':           [0.0027, False, 'NH'],
            'stretchDm31':    1000.0,
            'theta23':        [0.85, False],  
            "octant" :        "B",
            'theta13':        [0.155, True], 
            'mix_angle':      [0.95, True, 1.],
            'oscMode':        mode,
            'oscTables':      False,
            'baseline_llh':   0.0,
            # sterile neutrinos - fixing all, if want to fit - specify explicitly
            'dm41':           [1.0, True],
            'theta24':        [0.0, True],
            'theta34':        [0.0, True],
            # Normalizations and flux
            'norm':           [1., True], # Set to False only for your final result. Saves time.
            'norm_e':         [1., False],
            'norm_nc':        [1., True],
            'norm_tau':       [1., True], # Set to false for nutau apperance studies
            'gamma':          [0., False],  
            'nu_nubar':       [1., True],
            'nubar_ratio':     [0., False],
            'uphor_ratio':     [0., False],
            # Non-DIS cross sections
            'axm_qe':         [0., False],
            'axm_res':        [0., False],
            # Detector parameters
            'pid_bias':       [0., True],
            'hole_ice':       [0.02, False], 
            'hi_fwd':         [0., True],
            'domeff':         [1., False],
            'had_escale':     [1., True],
            # Backgrounds
            'atmmu_f':        [0.00, True, 'data'],
            'noise_f':        [0.0, True],
            # Other details
            'detector_syst':  True,
            'include_priors': True, 
            'printMode':      -1, # Select 0,1,2,3 for the iminuit print_levels. Select -1 for every step.
            "minimizer":     "Migrad",
            "llh_space":     "Poisson"
            }
## removing prior from gamma
priors = {"gamma": [0.01, 1000.0]}

fit_errors = {'error_dm31'   :     0.00001,
          'error_theta23':     0.0001,
          "error_atmmu_f":     0.0001,
          'error_norm_e':      0.0001,
          'error_norm':        0.0001,
          'error_gamma':       0.001,
          'error_domeff':      0.005, 
          'error_hole_ice':    0.001, 
          'error_axm_qe' :     0.001,
          'error_axm_res':     0.001,   
          'error_nubar_ratio':  0.001,
          'error_uphor_ratio':  0.001,   
          'error_theta24':     0.1 ,
          'error_theta34':     0.1 ,    
           }

import oscFit

fitter = oscFit.fitOscParams()


result = fitter(data_histograms = [tracks_data_histogram],
                data_loaders        = [tracks_loader],
                fit_settings        = fit_settings, 
                fit_priors          = priors,
                fit_errors          = fit_errors,
                store_fit_details   = False, 
                tol                 = 0.000001,
                evalOnly = False, 
                sin2fit = False)


f = open('onefit_output.pckl', 'w')
pickle.dump(result, f )
f.close()

