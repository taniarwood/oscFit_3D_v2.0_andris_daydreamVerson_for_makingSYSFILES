#################################################################
# oscfit_default_values.py
# Juan-Pablo Yanez (j.p.yanez@ualberta.ca)
# 1 September 2016
#
# This file defines multiple default values, both for the fitter
# and for the module that produces the histograms. Should be
# modified with care. Most options can be changed in the steering
# scripts.
#################################################################

from numpy import pi

# Parameters used in the fit
fit_keys = [
    # Standard oscillation parameters
    'dm31',
    'theta23',
    'theta13',
    # Standard oscillations in 2-flavor mode
    'mix_angle',
    # Sterile neutrino oscillations
    'dm41',
    'theta24',
    'theta34',
    # Normalizations of components
    'norm',

    # Tania normalizations
    'nu_frac1',
    'nu_frac2',
    'nu_frac3',
#    'nu_frac4',  QUESTION:  Is it hard coded that nu_frac4 has to be the one we tie to the normalizations.. shit. how hard-coded is this unfolding?
#suddenly i am thinking very.  I will need to look thorugh the code. 
    'nu_frac4',
    'nu_frac5',
    'nu_frac6',
    'nu_frac7',
    'nu_frac8', 

    'norm_tau',
    'atmmu_f',
    'noise_f',
    # Neutrino flux uncertainties
    'gamma',
    'norm_e',
    'nu_nubar',
    'nubar_ratio',
    'uphor_ratio',
    'nu_pi_scale',
    # Cross section uncertainties
    'norm_nc',
    'axm_qe',
    'axm_res',
    # Detection uncertainties
    'domeff',
    'hole_ice',
    'hi_fwd',
    'had_escale',
    'pid_bias'
    ]

# Data settings - input for data loader
default_data_settings = {
    'dm31':                     0.0025, 
    'theta23':                  0.74,
    'theta13':                  0.155,
    'mix_angle':                1., 
    'norm_nu':                  2.4., # In years
    'norm_e':                   1.,
    'norm_tau':                 1.,
    'norm_noise':               1.,
    'noise_f':                  None, #Set either a normalization or fraction
    'norm_nc':                  1.,
    'norm_atmmu':               0., 
    'atmmu_f':                  None, #Set either a normalization or fraction
    'nu_nubar':                 1., #
    'nubar_ratio':              0.,
    'uphor_ratio':              0.,
    'nu_pi_scale':              1.,
    'gamma':                    0.0, 
    'axm_qe':                   0.,
    'axm_res':                  0.,
    'dis_xa':                   0., # From Teppei's parameterization of DIS discrepancies
    'dis_xb':                   1., # Applied as b*x^(-a) ; x is the Bjorken-x
    'pid_bias':                 0.0,
    'domeff':                   1., 
    'hole_ice':                 0.02, 
    'hi_fwd':                   0.,
    'had_escale':               1.,
    'atmmu_template':           'data', #'mc',
    'simulation':               'baseline', 
    'oscMode':                  'TwoNeutrino',#'Vacuum',
    'oscTables':                False,
    'ma_variations':            True,
    'dis_reweight':             False, # Only turn it on if you have 'x' and 'scattering' in pickle files
    'add_detector_systematics': True,
    'correct_honda_flux':       False,
    'dm41':                     0.0,
    'theta24':                  0.0, 
    'theta34':                  0.0
    }

default_fit_settings  = {
    # Oscillations
    'simulation':     'baseline',
    'dm31':           [0.0025, False, 'NH'],
    'theta23':        [0.78, False], 
    'theta13':        [0.153, False],
    'mix_angle':      [0.95, False, 1.],
    'oscMode':        'TwoNeutrino',
    'oscTables':      False,
    'baseline_llh':   0.0,
    ### Sterile neutrinos - fixing all, if want to fit - specify explicitly
    'dm41':           [0.0, True],
    'theta24':        [0.0, True],
    'theta34':        [0.0, True],
    # Normalizations and flux
    'norm':           [1., True], # Set to False only for your final result. Saves time.
    'nu_frac1':         [0.1, False],
    'nu_frac2':         [0.1, False],
    'nu_frac3':         [0.1, False],
    'nu_frac4':         [0.1, False],
    'nu_frac5':         [0.1, False],
    'nu_frac6':         [0.1, False],
    'nu_frac7':         [0.1, False],
    'nu_frac8':         [0.1, False],

    'norm_e':         [1., False],
    'norm_nc':        [1., True],
    'norm_tau':       [1., True], # Set to false for nutau apperance studies
    'gamma':          [0., False],  
    'nu_nubar':       [1., False],
    'nubar_ratio':     [0., False],
    'uphor_ratio':     [0., True],
    'nu_pi_scale':     [1., False],
    # Non-DIS cross sections
    'axm_qe':         [0., True],
    'axm_res':        [0., True],
    # Detector parameters
    'pid_bias':       [0., True], # Assumes a linear PID bias. Think carefully before you use it.
    'hole_ice':       [0.02, False], 
    'hi_fwd':         [0., True],     # Set to False only if you know what you are doing.
    'domeff':         [1., False],
    'had_escale':     [1., True],
    # Backgrounds
    'atmmu_f':        [0.05, False, 'data'],
    'noise_f':        [0.0, True],
    # Other details
    'detector_syst':  True,
    'include_priors': True, 
    'printMode':      0, # Select 0,1,2,3 for the iminuit print_levels. Select -1 for every step.
    # Minimization
    'minimizer'     : "migrad", #used to be migrad 
#    'minimizer'      : "SLSQP",
    "octant"        : "B" , # B: both, L: left, R: right - A good fit should be run in both A and B.
    "llh_space"     : 'poisson',
    "blind_keys"    : [],
    "remove_blind_keys": True,# Sets the results from blind_keys to "BLIND", otherwise the results will be stored. 
    # False might be needed for studying difference between different years/sets, 
    # default is True - no keys from blind_keys will be printed and stored
    "fit_function"  : 'chi_squared', # Change to chi_squared for MSU-like behavior
            }

default_fit_priors  = {
    # Oscillations
    'theta13':        [0.155,0.008], # This error is on the sin(2*th13)**2 (PDG2014)
    # Normalizations and flux
    'norm_e':         [1., 0.05],
    'norm_nc':        [1., 0.2],
   # 'gamma':          [0.05, 0.1],
    'nubar_ratio':     [0., 1.],
    'uphor_ratio':     [0., 1.],
    'axm_qe':         [0., 1.],
    'axm_res':        [0., 1.],
    # Detector parameters
    'hole_ice':       [0.02, 0.01],
    #'hi_fwd':         [0., 10.],
    'domeff':         [1., 0.1],
    'had_escale':     [1., 0.1]
    }

# Expected errors of the fit - all need to be defined. Give a crude estimate
default_fit_errors = {        
    'dm31': 0.001,
    'theta23': 0.0001,
    'theta13': 0.008,
    'mix_angle': 0.1,
    'norm':  0.0001,
    'gamma': 0.1,

    'nu_frac1':0.4,
    'nu_frac2':0.4,
    'nu_frac3':0.4,
    'nu_frac4':0.4,
    'nu_frac5':0.4,
    'nu_frac6':0.4,
    'nu_frac7':0.4,
    'nu_frac8':0.2,

    'norm_e': 0.05,
    'norm_tau':0.1,
    'nu_nubar':0.6,
    'nubar_ratio':0.5,
    'uphor_ratio':0.5,
    'nu_pi_scale':0.1,
    'atmmu_f':0.01,
    'noise_f': 0.01,
    'axm_qe' :0.5,
    'axm_res':0.5,
    'norm_nc':0.1,
    'pid_bias':0.01,
    'domeff':0.02,
    'hole_ice':0.001,
    'hi_fwd':0.01,
    'had_escale':0.03,
    'theta24':0.2 ,
    'theta34':0.2 ,
    'dm41': 0.5
    }


default_fit_limits = {
    'theta23'        : (0., pi/2.), 
    'theta13'        : (0., pi/2.),
    'mix_angle'      : (0., 2.),
    'dm31'           : (0, 1.),
    'norm'           : (0.7, 1.3),
    'nu_frac1'         : (0., 1.),
    'nu_frac2'         : (0., 1.),
    'nu_frac3'         : (0., 1.),
    'nu_frac4'         : (0., 1.),
    'nu_frac5'         : (0., 1.),
    'nu_frac6'         : (0., 1.),
    'nu_frac7'         : (0., 1.),
    'nu_frac8'         : (0., 1.),


    'atmmu_f'        : (0.0, 1.0),
    'noise_f'        : (0.0, .2),
    'norm_e'         : (0.0, 4.0),
    'norm_tau'       : (0., 2.),
    'gamma'          : (-0.5, 0.5),
    'nu_nubar'       : (0, 1),
    'nubar_ratio'    : (-3, 3),
    'uphor_ratio'    : (-3, 3),
    'nu_pi_scale'    : (0, 3),
    'axm_qe'         : (-3, 3),
    'axm_res'        : (-3, 3),
    'domeff'         : (0.5, 1.5),
    'had_escale'     : (0.2, 1.2),
    'hole_ice'       : (0.00, 0.05),
    'hi_fwd'         : (-7., 7.),
    'norm_nc'        : (0.5, 1.5),
    'pid_bias'       : (-0.06, 0.06),
    'theta24'        : (0.00, pi/2.0), 
    'theta34'        : (0.00, pi/2.0),
    'dm41'           : (0, 10.),
    }
