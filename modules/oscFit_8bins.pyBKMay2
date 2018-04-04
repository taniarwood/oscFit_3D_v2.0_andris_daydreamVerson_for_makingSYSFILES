import numpy as np
from numpy.polynomial.polynomial import polyfit
import pickle
from copy import deepcopy
from sys import exit

from scipy.interpolate import interp1d
#import minuit2
import iminuit

import dataLoader, miscFunctions, oscfit_default_values
from likelihoods import likelihoods
reload(dataLoader)
reload(miscFunctions)

#np.set_printoptions(precision=0, linewidth=180)

class fitOscParams(object):

    def __init__(self):
        print '\n ************************************ '
        print ' ****** oscFit3D v1.1 - oscFit ******' 
        print ' ************************************\n'

        self.print_keys  = []

    def printMinimizationStep(self,LLH):
        print_list = ["%.5f" % LLH]
        for one_key in self.print_keys:
            if one_key == "dm31":
                print_list.append("%+.5f" % (self.hist_params[one_key] ))
            else: print_list.append("%+.5f" % self.hist_params[one_key])
        print_list[-1] = print_list[-1]
        print '  '.join(print_list)

    def printMinStepHeaders(self, fsettings):
        if self.printMode < 0:
            if len(self.print_keys) == 0:
                for key in oscfit_default_values.fit_keys:
                    if key in self.blind_keys: continue
                    if fsettings[key][1] == False: self.print_keys.append(key)
            print '\n' 
            print '-LLH   |'.rjust(14),
            for one_pKey in self.print_keys:
                print (one_pKey+" |").rjust(8),
        print ""

    # Functions needed for the scipy minimizers
    def constrain_All(self, p, lims):
        for tKey in lims:
            if tKey in p:
                tP, tLim = p[tKey], lims[tKey]
                if tP < tLim[0] or tP > tLim[1]: return False
            else:
                print "oscFit Warning: Got limit on parameter '", tKey,"' that does not exist in target function."
        return True

    def llhWrap(self, constrain_func, printMode, remArgs, **kwargs):
        if not constrain_func == None:
            is_ok = constrain_func(kwargs)
            if not is_ok: return is_ok
        tLlh = self.llhFunction(**kwargs)
        tInfo = "%0.5f"%(tLlh)+" "*(12-len(str(tLlh)))+" | "
        keys = kwargs.keys()
        keys.sort()
        for tK in keys:
            if tK in remArgs: continue
            if tK in self.blind_keys: continue
            tInfo += str(round(kwargs[tK],5)).rjust(9)
        if printMode: print tInfo
        return tLlh 

      
    # Likelihood function - I tried to make the arguments arbytrary (**kwargs). Minuit didn't like it.
    def llhFunction(self,dm31 = None,  
                         theta23 = None,  
                         theta13 = None, 
                         mix_angle = None,  
                         norm = None, 

                         nu_frac1 = None,
                         nu_frac2 = None,
                         nu_frac3 = None,
                         nu_frac4 = None,
			 nu_frac5 = None,
                         nu_frac6 = None,			 
			 nu_frac7 = None,
                         #nu_frac8 = None,

                         gamma = None,  
                         norm_e = None,  
                         norm_tau = None,  
                         nu_nubar   = None, 
                         nubar_ratio = None, 
                         uphor_ratio = None,
                         nu_pi_scale = None,
                         atmmu_f = None,  
                         noise_f = None, 
                         axm_qe = None,  
                         axm_res = None, 
                         norm_nc = None,  
                         pid_bias = None, 
                         domeff = None,  
                         had_escale = None,  
                         hole_ice = None, 
                         hi_fwd = None, 
                         dm41 = None,   
                         theta24 = None, 
                         theta34 = None):
        '''
        Returns the -LLH for a selection of values.
        '''

        # Passing the arguments to hist_params
        kwargs = locals()
        for iKey in oscfit_default_values.fit_keys:
            self.hist_params[iKey] = kwargs[iKey]

        # Obtaining the histograms from each loader
        neutrinos = []
        muons     = []
        noise     = []
        if self.fit_function == 'chi_squared':
            neutrinos_w2 = []
            muons_w2     = []
            noise_w2     = []
            
        nuperflavor = []  # Only used by the barlow likelihood
        barlow_out_singleh = [0]*len(self.data_loaders)
        barlow_out_totalh  = [0]*len(self.data_loaders)
        for index, one_loader in enumerate(self.data_loaders):
            nu_histo, mu_histo, noise_histo, perflavor  = one_loader.getSingleHistos(params = self.hist_params,
                                                                                     detailed_nu = True)
            neutrinos.append(nu_histo)
            muons.append(mu_histo)
            noise.append(noise_histo )
            nuperflavor.append(perflavor)    

            if self.fit_function == 'chi_squared':
                nu_histo_w2, mu_histo_w2, noise_histo_w2, perflavor_w2  =\
                                    one_loader.getSingleHistos(params = self.hist_params,
                                                               detailed_nu = True, weight_power = 2.)
                neutrinos_w2.append(nu_histo_w2)
                muons_w2.append(mu_histo_w2)
                noise_w2.append(noise_histo_w2)


        # Define the number of neutrinos expected and the global scaling factor to match observed events
        expected_events = self.data_events*norm
        nu_events       = expected_events*(1-atmmu_f-noise_f)
        norm_nu = nu_events/np.sum([np.sum(neutrinos[i]) for i in range(len(neutrinos))])

        # Split the total in the desired way
        norm_1 = nu_events*nu_frac1/(np.sum(neutrinos[0])*norm_nu)
        norm_2 = nu_events*nu_frac2/(np.sum(neutrinos[1])*norm_nu)
        norm_3 = nu_events*nu_frac3/(np.sum(neutrinos[2])*norm_nu)
        norm_4 = nu_events*nu_frac4/(np.sum(neutrinos[3])*norm_nu)
        norm_5 = nu_events*nu_frac5/(np.sum(neutrinos[4])*norm_nu)
        norm_6 = nu_events*nu_frac6/(np.sum(neutrinos[5])*norm_nu)
        norm_7 = nu_events*nu_frac7/(np.sum(neutrinos[6])*norm_nu)
#        norm_8 = nu_events*nu_frac8/(np.sum(neutrinos[7])*norm_nu)

#	norm_9 = nu_events*(1-nu_frac1-nu_frac2-nu_frac3-nu_frac4-nu_frac5-nu_frac6-nu_frac7-nu_frac8)/(np.sum(neutrinos[8])*norm_nu)
#        norm_9 = np.abs(norm_9)
        norm_8 = nu_events*(1-nu_frac1-nu_frac2-nu_frac3-nu_frac4-nu_frac5-nu_frac6-nu_frac7)/(np.sum(neutrinos[7])*norm_nu)
        norm_8 = np.abs(norm_8)


        # The minimizer does not know that norm_1,2,3 should sum to less than 1
        # Enforcing this by always using the absolute value of norm_4
        # The user should set more reasonable limits for the norms used by looking at MC prediction

        #nu_norms = [norm_1, norm_2, norm_3, norm_4,norm_5, norm_6, norm_7, norm_8, norm_9, 0.]
        nu_norms = [norm_1, norm_2, norm_3, norm_4, norm_5, norm_6, norm_7, norm_8, 0.]
        #print 'Expected events', expected_events, norm_nu
        norm_atmmu = 0.
        norm_noise = 0.
        if atmmu_f != 0:
            atmmu_sum = np.sum([np.sum(x) for x in muons])
            if atmmu_sum > 0:
                norm_atmmu = expected_events*atmmu_f/atmmu_sum
        if noise_f != 0:
            noise_sum = np.sum([np.sum(x) for x in noise])
            if noise_sum > 0:
                norm_noise = expected_events*noise_f/noise_sum
        #print 'Atmmu_norm', atmmu_sum, norm_atmmu

        # Producing the expectation for each loader using normalizations. Adding up the LLH contributions
        # Even if you use chi-squared, I'm still calling the output LLH for the sake of not renaming stuff
        LLH = 0
        ref_histo = np.zeros_like(self.exp_histograms[0])
        ref_histo_w2 = np.zeros_like(ref_histo)

        for index, loader in enumerate(self.data_loaders):
            ref_histo += neutrinos[index]*norm_nu*nu_norms[index] +\
                         muons[index]*norm_atmmu + noise[index]*norm_noise
            print 'Nu norms index', nu_norms[index],np.sum( neutrinos[index]*norm_nu*nu_norms[index])

            if self.fit_function == 'chi_squared':
                ref_histo_w2 += neutrinos_w2[index]*norm_nu**2*nu_norms[index]**2 +\
                                muons_w2[index]*norm_atmmu**2 + noise_w2[index]*norm_noise**2 +\
                                (loader.atmmu_histo['data'] + 
                                 (loader.atmmu_histo['data'] - loader.atmmu_histo['data_aux'])**2)*norm_atmmu**2

        # Reference histogram
        print ref_histo
        print np.sum(ref_histo)

        nonzero           = ref_histo > 0

        if self.fit_function == 'llh':
            if self.llh_space == "poisson":
                LLH += sum(ref_histo[nonzero] - self.exp_histograms[0][nonzero]*np.log(ref_histo[nonzero]))
		#makes liklihood values easier to read. now '0' is a 'perfect' fit. 
		#can go to chi2 more easily 
		#this is a likeli of MC fitting data, term added is llh that data fits itself. 
		#makes this a delta LLh so wilks etc works.
                
                # nonzero = self.exp_histograms[0] > 0
                # LLH -= sum(self.exp_histograms[0][nonzero] - 
                #            self.exp_histograms[nonzero]*np.log(self.exp_histograms[nonzero]))

            elif self.llh_space == "barlow":
                print 'Barlow and Tania do not get along yet'
                sys.exit()
                # Use the likelihood class from Likelihood.py to compare the histograms
                # Histograms merged into one structure for simplicity
                keys = ['numu_histo','nue_histo','nc_histo','nutau_histo']
                histograms = [nuperflavor[index][key] * norm_nu for key in keys]
                unweighted_histograms = [loader.unweighted_histograms[key]*self.barlow_scaling['nu'] 
                                         for key in keys]
                
                # In some cases, the atmospheric muon rate may drop to 0 when the fitter hits norm_atmmu=0. 
                # We should protect against that case here.
                # Scaling raw number of events. Avoiding having empty bins
                if np.sum(muons[index] * norm_atmmu):
                    histograms.append(muons[index] * norm_atmmu)
                    if self.barlow_scaling['muon'] != 1:
                        zerobins = loader.unweighted_histograms['atmmu_histo'] == 0
                        loader.unweighted_histograms['atmmu_histo'][zerobins] = 0.5
                    unweighted_histograms.append(loader.unweighted_histograms['atmmu_histo'] * 
                                                 self.barlow_scaling['muon'])

                if np.sum(noise[index]):
                    histograms.append(noise[index] * norm_noise)
                    if self.barlow_scaling['noise'] != 1:
                        zerobins = loader.unweighted_histograms['noise_histo'] == 0
                        loader.unweighted_histograms['noise_histo'][zerobins] = 0.5
                    unweighted_histograms.append(loader.unweighted_histograms['noise_histo'] *
                                                 self.barlow_scaling['noise'])


                # Numpyify them
                histograms = np.array(histograms)
                unweighted_histograms = np.array(unweighted_histograms)

                # The Barlow likelihood REQUIRES that every bin have an estimation of the weight/event. If we
                # want to use the Barlow method, set the weight/event in the bins with 0 events equal to the
                # minimum nonzero weight/event of that histogram
                histograms /= unweighted_histograms
                histograms = np.nan_to_num(histograms)
                for i in range(histograms.shape[0]):
                    indicies = histograms[i]==0
                    if not np.sum(indicies) == 0:
                        m = histograms[i][ np.logical_not(indicies) ]

                        # Weird case: every bin has a total weight of 0. This can happen if you have 0
                        # of one of the neutrino types [eg, taus]
                        if np.sum(m) == 0: continue
                        histograms[i][indicies] = np.min( m )

                # Set up the likelihood object with the histograms you want to compare
                likelihood_obj = likelihoods()
                likelihood_obj.SetData( self.exp_histograms[index] )
                likelihood_obj.SetMC( histograms )
                likelihood_obj.SetUnweighted( unweighted_histograms )

                # Run it with the requested LLH space.
                LLH += likelihood_obj.GetLLH( self.llh_space )

                # Barlow histograms cannot be retrieved later. Need to be kept in memory.
                barlow_out_singleh[index] = likelihood_obj.GetSinglePlots()
                barlow_out_totalh[index]  = likelihood_obj.GetPlot()

            else: 
                print "Unknown LLH space"
                LLH = np.nan       

        elif self.fit_function == 'chi_squared':
            LLH += sum( (ref_histo[nonzero] - self.exp_histograms[0][nonzero])**2/
                     (self.exp_histograms[0][nonzero] + ref_histo_w2[nonzero]))

        # Adding the priors after the bin-wise LLH contributions
        if self.include_priors:
            for param in self.fit_priors:
                # Learn about LLH vs chi-squared, y'all!
                if self.fit_function == 'chi_squared':
                    prior_factor = 1.
                elif self.fit_function == 'llh':
                    prior_factor = 0.5

                
                LLH += prior_factor*((self.hist_params[param]-self.fit_priors[param][0])/\
                                     self.fit_priors[param][1])**2


        # Store the output histograms of Barlow LLH - cannot be recovered later
        if self.llh_space == 'barlow' and LLH < self.barlow_bestfit_llh:
            self.barlow_bestfit_llh = LLH
            self.barlow_histograms     = barlow_out_singleh
            self.barlow_histograms_tot = barlow_out_totalh

        if np.isnan(LLH):
            print 'oscFit: The histogram requested returned NAN values. Check!'
            print ref_histo
            print self.hist_params
            exit()

        if self.printMode < 0:               
            self.printMinimizationStep(LLH)

        return LLH
             

    def __call__(self, data_histograms = [], 
                 data_loaders         = [], 
                 fit_settings         = {},
                 fit_priors           = {},
                 fit_errors           = {}, 
                 fit_limits           = {},
                 ncalls               = 1000, # For migrad. Simplex will use twice the number.
                 store_fit_details    = False,
                 return_minuit_object = False,
                 tol                  = -1.0,
                 evalOnly             = False,
                 barlow_scaling       = {'nu':1., 'muon':1., 'noise':1.}
                 ):
        '''
        Runs the fit and returns the values at the best fit point.
        Set store_fit_details to False for likelihood scans. Set it to True to obtain covariance, correlations.

        Parameters
        ----------

        data_histograms - Binned data
        
        data_loaders - dataLoader instances to obtain MC histograms

        fit_settings - Dictionary containing the initial values, fixed parameters and other settings

        fit_priors - Dictionary containing the priors to be used. Not to be toyed around with.
        
        fit_errors - Dictionary containing the errors (estimated distance to minimum) for parameters,
                     in format: error_parname, example - {'error_dm31':0.0001 }
        
        store_fit_details - Store covariance and correlation matrices. Disable for scans.

        return_minuit_object - Return the iminuit. Use in interactive sessions. The object cannot be stored.
        
        tol  - set tolerance to the minimizer, if < 0.0 - default will be used (0.01 so far)

        evalOnly - only evaluate LLH without doing fit
        '''

        # Verifying the input makes sense
        if len(data_loaders) == 0:
            print 'oscFit: At least one dataLoader object has to be given as an argument'
            exit()
        miscFunctions.checkDictParams(oscfit_default_values.default_fit_settings,
                                      fit_settings, 'oscFit')
        fsettings = deepcopy(oscfit_default_values.default_fit_settings)
        fsettings.update(fit_settings)
        ferrors = deepcopy(oscfit_default_values.default_fit_errors)
        ferrors.update(fit_errors)
        miscFunctions.checkDictParams(oscfit_default_values.default_fit_limits,
                                      fit_limits, 'oscFit')
        flimits = deepcopy(oscfit_default_values.default_fit_limits)
        miscFunctions.checkDictParams(oscfit_default_values.default_fit_settings, 
                                      fit_priors, 'oscFit')
        self.fit_priors = deepcopy(oscfit_default_values.default_fit_priors)
        self.fit_priors.update(fit_priors)

        # These parameters do not change after every step and are needed in dataLoader
        self.hist_params = {'atmmu_template'            : fsettings['atmmu_f'][2],
                            'add_detector_systematics'  : fsettings['detector_syst'],
                            'simulation'                : fsettings['simulation'],
                            'oscMode'                   : fsettings['oscMode'],
                            'oscTables'                 : fsettings['oscTables']}

        # Gathering some additional info
        print 'oscFit: Fitter running in', fsettings['oscMode'], 'mode'
        self.exp_histograms = data_histograms
        self.data_events    = np.sum([np.sum(x) for x in self.exp_histograms])
        print 'oscFit:', self.data_events, ' events in data'

        self.include_priors = fsettings['include_priors']
        self.printMode      = fsettings['printMode']
        self.blind_keys     = fsettings['blind_keys']
        self.data_loaders   = data_loaders
        self.llh_space      = fsettings['llh_space'].lower()
        self.curFitOctant   = fsettings["octant"] 
        self.remove_blind_keys    = fsettings['remove_blind_keys']

        self.fit_function  = fsettings['fit_function']

        # Declare variables needed by the Barlow LLH
        if self.llh_space == 'barlow':
            self.barlow_scaling = barlow_scaling
            self.barlow_bestfit_llh = np.inf

        # If there are blind keys, need some protection
        if len(self.blind_keys)>0: 
            print "oscFit: the following keys will be blind:"
            print "\t",self.blind_keys
            print "Setting dataLoader.verbose=False"
            for one_loader in self.data_loaders:
                one_loader.verbose=False

        # Define the hierarchy, thus the sign of dm31
        if fsettings['dm31'][2] == 'NH':
            print 'oscFit: Fitting in normal hierarchy mode'
            flimits['dm31'] = (0, 7E-3)
            fsettings['dm31'][0] = abs(fsettings['dm31'][0])
        elif fsettings['dm31'][2] == 'IH':
            print 'oscFit: Fitting in inverted hierarchy mode'
            flimits['dm31'] = (-7E-3, 0)
            fsettings['dm31'][0] = -1.*abs(fsettings['dm31'][0])
        else:
            print 'oscFit: Fit settings for dm31 are not valid - ', fsettings['dm31']
            exit()

        # Define the octant (it can be an issue)
        flimits['theta23'] = ( 0.0       , np.pi/2.0 )
        if "octant" in fsettings:
            print "oscFit: found octant in settings",fsettings["octant"]
            if fsettings["octant"] == "L":      flimits['theta23'] = ( 0.0       , np.pi/4.0 )
            elif fsettings["octant"] == "R":    flimits['theta23'] = ( np.pi/4.0 , np.pi/2.0 )
            else:                               flimits['theta23'] = ( 0.0       , np.pi/2.0 )
            print "found octant constrain to theta23 as "+str(flimits['theta23'])

        # Define the mixing angle limit
        flimits['mix_angle'] = (0.,fsettings['mix_angle'][2])

        # Evaulate the LLH once with the input values
        starting_eval_dict = {}
        for one_key in oscfit_default_values.fit_keys:
            starting_eval_dict[one_key] = fsettings[one_key][0]                                
	print starting_eval_dict
        llh_testval = self.llhFunction(**starting_eval_dict)

        print 'oscFit: Testing the LLH function ', llh_testval
        if evalOnly: 
            print "Only evaluation of the LLH was selected, returnung LLH value"
            return { "llh" : llh_testval } 

        self.print_keys = []
        
        print 'oscFit: Fit settings'
        for one_key in fsettings:
            print '\t', one_key, '\t',fsettings[one_key]
        print ''


        # Passing initial values. Defining what gets fixed.
        if not self.hist_params['add_detector_systematics']:
            fsettings['domeff'][1]     = True
            fsettings['hole_ice'][1]   = True
            fsettings['hi_fwd'][1]     = True
            fsettings['had_escale'][1] = True
        # Keeping the possibility to fit in 2-flavor mode (for checks)
        if self.hist_params['oscMode'] == 'TwoNeutrino':
            fsettings['theta13'][1]   = True
            fsettings['theta23'][1]   = True
        else:
            fsettings['mix_angle'][1] = True


        #Running minimizer
        if fsettings["minimizer"].lower() == "migrad":        

            # Building the arguments
            kwargs = {}

            # Passing the fit PARAMETERS, ERRORS, LIMITS and FIX settings to kwargs
            for one_var in oscfit_default_values.fit_keys:
                kwargs[one_var] = fsettings[one_var][0]
                kwargs['error_'+one_var] = ferrors[one_var]
                kwargs['fix_'+one_var] = fsettings[one_var][1]
                kwargs['limit_'+one_var] = flimits[one_var]

            kwargs['print_level'] = self.printMode

            if fsettings['fit_function'] == 'chi_squared':
                kwargs['errordef']    = 1.
            elif fsettings['fit_function'] == 'llh':
                kwargs['errordef']    = 0.5

            m = iminuit.Minuit(self.llhFunction, **kwargs)

     
            m.strategy  = 1
            m.errordef  = 0.5
            if tol <= 0.0: tol = 0.01 # Zero also does not make sence
            print "oscFitL: Tolerance for iminuit: ", tol
            m.tol       = tol

            # Starting!
            good_fit = False
            # Let's try to fit at least 3 times. If the fit does not converge at all, go with simplex.
            fit_counter = 0
            routine     = 'migrad'
            while (not good_fit) and (fit_counter < 3):
                fit_counter += 1
                try:
                    print 'oscFit: Firing up the minimizer!'
                    m = iminuit.Minuit(self.llhFunction, **kwargs)
                    m.strategy  = 1
                    m.errordef  = 0.5
                    m.tol       = tol
                    # Minimize
                    self.printMinStepHeaders(fsettings)
                    m.migrad(ncall=ncalls)
                    good_fit = True
                except:
                    # If the minimizer didn't converge, modify the errors
                    print 'oscFit: MIGRAD did not converge. Trying again with different initial step-sizes'
                    for kwargs_key in kwargs:
                        if 'error_' in kwargs_key:
                            kwargs[kwargs_key] *= 10**np.random.normal()

            if not good_fit:
                try:
                    routine = 'simplex'
                    self.printMinStepHeaders(fsettings)
                    m = iminuit.Simplex(ncall=ncalls*2)
                    good_fit = True
                except:
                    print 'oscFit: SIMPLEX failed as well'

            if good_fit:
                if store_fit_details:
                    try:
                        m.hesse()
                        hesse_errors = True
                    except:
                        print 'oscFit: Tried hesse errors and failed'
                        hesse_errors = False
            else:
                print "Failed migrad with:", m.values

            # Storing the results
            results = {'llh':m.fval, 
                       'expected_events':[np.sum(x) for x in self.exp_histograms], 
                       'atmmu_template':self.hist_params['atmmu_template'], 
                       'oscMode':self.hist_params['oscMode'], 
                       'fit_function':self.fit_function,
                       'add_detector_systematics':self.hist_params['add_detector_systematics']}
	    
            if good_fit and store_fit_details and hesse_errors:
                results['covariance']  = m.covariance
                results['parameters']  = m.parameters
                results['errors']      = m.errors

                try:
                    results['corr_matrix'] = np.array(m.matrix(correlation=True))
                    results['hesse_errors'] = hesse_errors
                except:
                    print 'oscFit: HESSE appears to have failed. That is ok! '
                    print 'You still have the fit and we do not use the errors anyway.'
                    print 'This is telling you that \n1) the error estimate will not be good and'
                    print '2) No correlation matrix will be available. You can live with that'
                    results['hesse_errors'] = False

            results.update(m.values)

        else:
            print "oscFit: scipy minimizer is selected"
            try:     
                from scipy import optimize 
            except ImportError: 
                print "scipy.optimize could not be imported. Stop minimization here." 
                exit() 
            import minimizerSettings
            known_minimizers  = minimizerSettings.known_minimizers
            minimizer_options = minimizerSettings.minimizer_options
            if not fsettings["minimizer"] in known_minimizers: 
                print "oscFit: unknown minimizer", fsettings["minimizer"]
                exit()

            # convert stuff such that it is suitable for the scipy.optimize syntax # 
            tArgs   = deepcopy(oscfit_default_values.fit_keys)
            tDirec  = ferrors
            tLimits = flimits

            fixed = {}
            
            # remove unnecessary args and convert limits and directions #
            remArgs = []
            print "Arguments used: "
            for tArg in tArgs:
                print tArg,"   ", fsettings[tArg][1]
                if fsettings[tArg][1]:
                    remArgs.append(tArg)
                    del flimits[tArg]
                    del ferrors[tArg]
                    fixed.update(  { tArg : fsettings[tArg][0] } )
            for tArg in remArgs: tArgs.remove(tArg)

            # JP wtf is there a function defined in here? - move it to the minimizers file
            # constrain fct included in llh fct for non-bounded optimizers to operate on bounded spaces #
            # def constrain_All(p, lims):
            #     for tKey in lims:
            #         if tKey in p:
            #             tP, tLim = p[tKey], lims[tKey]
            #             if tP < tLim[0] or tP > tLim[1]: return False
            #         else:
            #             print "Warning: Got limit on parameter '", tKey,"' that does not exist in target function."
            #     return True
     
            # JP again, why is this even here!?
            # wrap llh function, get more ifno during minimization & allow boundaries for unbounded minimizers #
            # def llhWrap(constrain_func, printMode, **kwargs):
            #     if not constrain_func == None:
            #         is_ok = constrain_func(kwargs)
            #         if not is_ok: return is_ok
            #     tLlh = self.llhFunction(**kwargs)
            #     tInfo = "%0.5f"%(tLlh)+" "*(12-len(str(tLlh)))+" | "
            #     keys = kwargs.keys()
            #     keys.sort()
            #     for tK in keys:
            #         if tK in remArgs: continue
            #         if tK in self.blind_keys: continue
            #         tInfo += str(round(kwargs[tK],5)).rjust(9)
            #     if printMode: print tInfo
            #     return tLlh

            # prepare fct, seed and bounds #
            tFunc       = lambda p: self.llhWrap( constrain_func, fsettings['printMode'], remArgs, 
                                             **{ tA:tp for tA, tp in zip(tArgs, p) + fixed.items() } )
            tX0         = [  fsettings[tA][0] for tA in tArgs  ]
            tBounds     = [  tLimits[tA]      for tA in tArgs  ]
            tDirec      = [  np.roll( [tDirec[tA] ] + [0.0]*(len(tDirec)-1) , r) for r, tA in enumerate(tArgs)  ]
            
            if fsettings["minimizer"] in ("Nelder-Mead", "Powell"):    
                constrain_func = lambda p: self.constrain_All(p, tLimits) 
                if fsettings["minimizer"] == "Powell":
                    minimizer_options[fsettings["minimizer"]]["options"].update(  { "direc" : tDirec }  )
            else: constrain_func = None

            print "Found - Args:",len(tArgs),"Bounds:",len(tBounds)," Direcs:",len(tDirec),"Limits:",len(tLimits)
            
            # print header #
            def print_header():
                tHeader     = "LLH:"+" "*(12-len("LLH:"))+" | "
                allKeys = tArgs+fixed.keys()
                allKeys.sort()
                for tA in allKeys: 
                    if tA in self.blind_keys: continue
                    if tA in remArgs: continue
                    tHeader += str(tA)[0:7].rjust(9)
                print tHeader
            if fsettings['printMode']: print_header()

            # actual minimization routine #
            good_fit    = True
            #try:
            if fsettings["minimizer"] in ("Nelder-Mead", "Powell"):
               res = optimize.minimize(tFunc, tX0, **minimizer_options[fsettings["minimizer"]])
            else:
               res = optimize.minimize(tFunc, tX0, bounds=tBounds, **minimizer_options[fsettings["minimizer"]])
            if not res["success"]: good_fit = False
            
            # ... as a reminder #
            if fsettings['printMode']: print_header()
            
            # Storing the results convert results #
            results = {'llh': res["fun"], 
                       'expected_events':[np.sum(x) for x in self.exp_histograms], 
                       'atmmu_template':self.hist_params['atmmu_template'], 
                       'oscMode':self.hist_params['oscMode'], 
                       'add_detector_systematics':self.hist_params['add_detector_systematics']}
                       
            if good_fit and store_fit_details:
                results['covariance']  = None
                results['parameters']  = tArgs
                results['corr_matrix'] = None
                
            results.update(  { tKey : res["x"][i] for i, tKey in enumerate(tArgs) }  )
            results.update(  fixed  )

            routine = fsettings["minimizer"]
            # as before - no different treatment due to the minimizers any more # 
            
        results['chi_squared']        = 0
        results['expected_events_mc'] = []
        results['simulation']         = self.hist_params['simulation']
        results['atmmu_template']     = self.hist_params['atmmu_template']
        results['oscMode']            = self.hist_params['oscMode']
        results['add_detector_systematics'] = self.hist_params['add_detector_systematics']
        results['successful_fit']     = good_fit
        results['fit_settings']       = fsettings
        results['min_routine']        = routine
        results['dm31'] = results['dm31']
        if results['successful_fit'] and results.has_key('errors') and ('dm31' in results['errors'].keys()):
            results['errors']['dm31'] =results['errors']['dm31']

        print '\n\n ************ FIT FINISHED ************\n'
        # ... oook, except here:
        if fsettings["minimizer"].lower() == "migrad":
            print 'Migrad EDM\t', m.edm
            print 'Migrad NCalls\t', m.ncalls
        else:
            print 'Calls to ', fsettings["minimizer"], " target fct.: \t", res["nfev"]

        # Fit is over. Retrieve the parameters and calculate normalizations (based on fractions fit)

        neutrinos = []
        muons     = []
        noise     = []
        for index, one_loader in enumerate(self.data_loaders):
            nu_histo, mu_histo, noise_histo  = one_loader.getSingleHistos(results)
            neutrinos.append(nu_histo)
            muons.append(mu_histo)
            noise.append(noise_histo)

        # Obtaining the normalizations that correspond to the fractions fit before
        expected_events_fit = self.data_events*results['norm']
        nu_events = expected_events_fit*(1-results['atmmu_f']-results['noise_f'])
        results['norm_nu'] = nu_events/(np.sum([np.sum(neutrinos[i])for i in range(len(neutrinos))])*\
                                        self.data_loaders[0].sec2years)
        norm_nu_secs = results['norm_nu']*self.data_loaders[0].sec2years
        nu_norms  = [nu_events*results['nu_frac1']/(np.sum(neutrinos[0])*norm_nu_secs),
                     nu_events*results['nu_frac2']/(np.sum(neutrinos[1])*norm_nu_secs),
                     nu_events*results['nu_frac3']/(np.sum(neutrinos[2])*norm_nu_secs),
		     nu_events*results['nu_frac4']/(np.sum(neutrinos[3])*norm_nu_secs),
                     nu_events*results['nu_frac5']/(np.sum(neutrinos[4])*norm_nu_secs),
                     nu_events*results['nu_frac6']/(np.sum(neutrinos[5])*norm_nu_secs),
                     nu_events*results['nu_frac7']/(np.sum(neutrinos[6])*norm_nu_secs),
        #             nu_events*results['nu_frac8']/(np.sum(neutrinos[7])*norm_nu_secs),
                     nu_events*(1-results['nu_frac1']-results['nu_frac2']-results['nu_frac3']-results['nu_frac4']-results['nu_frac5']-results['nu_frac6']-results['nu_frac7'])/\
                     (np.sum(neutrinos[7])*norm_nu_secs)]

        results['nu_norms']   = nu_norms
        results['norm_atmmu'] = 0.
        if results['atmmu_f'] != 0:
            results['norm_atmmu'] = expected_events_fit*results['atmmu_f']/np.sum([np.sum(x) for x in muons])
            if self.hist_params['atmmu_template'] != 'data':
                results['norm_atmmu'] /= self.data_loaders[0].sec2years
        results['norm_noise'] = 0
        if results['noise_f'] != 0:
            results['norm_noise'] = expected_events_fit*results['noise_f']/np.sum([np.sum(x) for x in noise])  
            

        # All of this changes for Tania
        final_ref_histo = np.zeros_like(self.exp_histograms[0])
        nonzero_bins = 0
        for index, one_loader in enumerate(self.data_loaders):
            # Barlow LLH has the parameters stored. Poisson LLH does not.
            if self.llh_space == 'poisson':
                final_histogram   = one_loader.getRefHisto(results)
            # Assume that you pass the bkg loader last always
            if index < len(self.data_loaders)-1:
                final_histogram *= nu_norms[index]

            final_ref_histo += final_histogram


            # elif self.llh_space == 'barlow':
            #     final_histogram   = self.barlow_histograms_tot[index]

        results['expected_events_mc'].append(np.sum(final_ref_histo))

        # Calculate the chi-squared - MSU modifications would go in here
        nonzero = (self.exp_histograms[0] > 0)*(final_ref_histo > 0) 
        results['chi_squared'] += np.sum((final_ref_histo[nonzero]-
                                          self.exp_histograms[0][nonzero])**2/
                                         self.exp_histograms[0][nonzero])

        nonzero_bins += np.sum(nonzero)
        print '\noscFit: Loader ', index, ' nonzero bins: ', np.sum(nonzero)

        if self.llh_space == 'barlow':
            results['barlow_histograms']     = self.barlow_histograms
            results['barlow_histograms_tot'] = self.barlow_histograms_tot

        results['nonzero_bins'] = nonzero_bins
        print 'oscFit: Total nonzero bins: ', np.sum(nonzero_bins)
        
        # Print the result (don't print the matrices!)
        print '\n **** oscFit3D results ****'
        skip_keys = ['fit_settings','corr_matrix','covariance', 'parameters', 
                     'errors', 'barlow_histograms', 'barlow_histograms_tot']
        result_keys = results.keys()
        result_keys.sort()
        for one_key in result_keys:
            if one_key in skip_keys:
                continue
            elif one_key in self.blind_keys: 
                if self.remove_blind_keys:
                    results[one_key]= "BLIND"
                print one_key, '\t BLIND KEY'
                continue
            if results['successful_fit'] and results.has_key('errors') and (one_key in results['errors'].keys()):
                print one_key, '\t',  results[one_key], '\t', results['errors'][one_key]
            else:
                print one_key, '\t',  results[one_key]
        print '\n'
        if good_fit:
            print 'The fit was successful!\n'
        else:
            print 'The fit was NOT successful. Check the steps to figure out why.\n'

        if return_minuit_object:
            return results, m
            # You can access the minuit object when running the fit interactively (ipython)
        else:
            return results
