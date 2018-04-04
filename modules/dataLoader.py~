"""
.. module:: dataLoader
   :synopsis: Machinery to load all data files into memory and produce histograms.

.. moduleauthor:: Juan-Pablo Yanez (j.p.yanez@ualberta.ca)


"""


import os, sys, pickle
settings_dir = os.path.dirname(os.path.realpath(__file__)).split('/modules')[0] + '/user_settings'
sys.path.append(settings_dir)
import numpy as np
from copy import deepcopy


import oscCalc
reload(oscCalc)

from oscCalc import propagationDistance
from miscFunctions import mergeDictionaries, addStatisticalFluctuations
import miscFunctions as mf
import systematicFunctions as sf
import oscfit_default_values
reload(oscfit_default_values)

# The reloading is necessary for testing, interactive use.
reload(sf)
reload(mf)
np.set_printoptions(precision=3, linewidth=120)


class dataLoader(object):
    
    def addCuts(self, in_dict, keep_observables = False, quantile_mc = (0., 1.)):
        '''
        Adds cuts while loading data & MC. Events outside the analysis range (reco energy & zenith) 
        are removed by default, unless you activate the keep_observables flag.
        Additional cuts can be placed for variables in the pickle file (e.g. vertex_z).
        The additional cuts are defined in the "extra_cuts" variable of dataLoader.
        '''

        these_cuts = deepcopy(self.extra_cuts)
        if keep_observables:
            for one_obs in self.observables:
                these_cuts.pop(one_obs)

        # If the data is loaded in L/E the variable has to be defined
        if 'LE' in self.observables:
            in_dict['LE'] = np.log10(propagationDistance(in_dict['reco_zenith'])/in_dict['reco_energy'])

        sel_bool = np.ones(len(in_dict[self.observables[-1]]), dtype=bool)

        # The quantile is applied on the sel_bool by setting it to False
        mi = int(quantile_mc[0]*len(sel_bool))
        ma = int(quantile_mc[1]*len(sel_bool))
        sel_bool[:mi] = False
        sel_bool[ma:] = False

        for cut_key in these_cuts.keys():
            if not in_dict.has_key(cut_key):
                print 'dataLoader: ERROR! Applying a cut on a non-existant variable ', cut_key
                sys.exit()
            sel_bool *= (in_dict[cut_key] >= np.min(these_cuts[cut_key])) *\
                        (in_dict[cut_key] < np.max(these_cuts[cut_key]))

        new_dict = {}

        # Producing the new dictionaries (nested CC and NC for MC)
        if 'interaction' in in_dict.keys():
            new_dict['CC'] = {}
            new_dict['NC'] = {}
            ccbool = in_dict['interaction'] == 1
            ncbool = in_dict['interaction'] == 2
            for dict_keys in in_dict.keys():
                new_dict['CC'][dict_keys] = in_dict[dict_keys][sel_bool*ccbool]
                new_dict['NC'][dict_keys] = in_dict[dict_keys][sel_bool*ncbool]
        else:
            for dict_keys in in_dict.keys():
                new_dict[dict_keys] = in_dict[dict_keys][sel_bool]

        return new_dict

        
    def countMCevents(self):
        return (len(self.nue['baseline']['CC'][self.hist_obs[0]])+
                len(self.nue['baseline']['NC'][self.hist_obs[0]])+
                len(self.numu['baseline']['CC'][self.hist_obs[0]])+
                len(self.numu['baseline']['NC'][self.hist_obs[0]])+
                len(self.nutau['baseline']['CC'][self.hist_obs[0]]))

    def loadAtmMuMC(self, baseline_only = False, use_kde = False, quantile_mc = (0., 1.)):
        '''
        Loads atmospheric muon mc. Similar to the neutrino function, but simpler
        '''

        systKeys = self.user.atmmu_sets.keys()
        if len(systKeys) == 0:
            return
        if baseline_only:
            systKeys = ['baseline']


        # Loading MC sets
        for set_key in systKeys:
            one_set = self.user.atmmu_sets[set_key]
            for set_subkey in one_set.keys():
                one_subset = one_set[set_subkey]
                if self.atmmu.has_key(one_subset): continue
                self.atmmu[one_subset] = self.addCuts(pickle.load(open(os.path.join(
                    self.user.data_dir, one_subset + '.pckl'))), use_kde, quantile_mc)
        self.atmmu['baseline'] = self.atmmu[self.user.atmmu_sets['baseline']['baseline']]
        # JPmod Table preparation would go in here
        

    def loadNeutrinoMC(self, baseline_only = False, use_kde = False, quantile_mc = (0.,1.)):
        '''
        Loads neutrino MC sets and applies the cuts defined. GENIE and NuGen are mixed at the break_energy 
        to create a "fixed_baseline" set. NuGen MC is scaled to match GENIE. Only necessary for old GENIE MC.
        '''
        if baseline_only:
            systKeys = ['baseline']
        else:
            systKeys = self.user.mc_sets.keys()
        
        # Loading MC sets
        for set_key in systKeys:
            one_set = self.user.mc_sets[set_key]
            for set_subkey in one_set.keys():
                one_subset = one_set[set_subkey]

                # Check if the MC set has been already loaded.
                if self.nue.has_key(one_subset): continue
                self.nue[one_subset]   = self.addCuts(pickle.load(open(self.user.genie_p1['nue']+
                                                                       one_subset+self.user.genie_p2)), use_kde,
                                                      quantile_mc)
                self.numu[one_subset]  = self.addCuts(pickle.load(open(self.user.genie_p1['numu']+
                                                                       one_subset+self.user.genie_p2)), use_kde,
                                                      quantile_mc)
                self.nutau[one_subset] = self.addCuts(pickle.load(open(self.user.genie_p1['nutau']+
                                                                       one_subset+self.user.genie_p2)), use_kde,
                                                      quantile_mc)


        self.nue['baseline']   = self.nue[self.user.mc_sets['baseline']['baseline']]
        self.numu['baseline']  = self.numu[self.user.mc_sets['baseline']['baseline']]
        self.nutau['baseline'] = self.nutau[self.user.mc_sets['baseline']['baseline']]
        # without this it would be hard to keep track off where the baseline came from # 
        self.baseline_key      = self.user.mc_sets['baseline']['baseline'] 

        # Producing the "fixed_baseline" in case there are nugen sets available
        if hasattr(self.user, 'nugen_nue') and hasattr(self.user, 'nugen_numu') and self.break_energy > 0:
            if self.user.nugen_nue != None and self.user.nugen_numu != None:
                # Loading NUGEN sets to fix the baselines (HE neutrinos missing in GENIE)
                nue_nugen = self.addCuts(pickle.load(open(self.user.nugen_nue)), use_kde, quantile_mc)
                numu_nugen = self.addCuts(pickle.load(open(self.user.nugen_numu)), use_kde, quantile_mc)
                self.nue['fixed_baseline']  = mf.mergeDictionaries(genie_dict = self.nue['baseline'], 
                                                                   nugen_dict = nue_nugen, 
                                                                   nugen_ccfactor = self.user.nugen_nuecc,
                                                                   nugen_ncfactor = self.user.nugen_nuenc,
                                                                   weight_keys = [self.we, self.wmu, self.wmu_pi],
                                                                   break_energy = self.break_energy)

                self.numu['fixed_baseline'] = mf.mergeDictionaries(genie_dict = self.numu['baseline'], 
                                                                   nugen_dict = numu_nugen,
                                                                   nugen_ccfactor = self.user.nugen_numucc,
                                                                   nugen_ncfactor = self.user.nugen_numunc,
                                                                   weight_keys = [self.we, self.wmu, self.wmu_pi],
                                                                   break_energy = self.break_energy)


                self.nutau['fixed_baseline'] = self.nutau['baseline']

        # JPmod table preparation should include also NC (if I move the flux to tables as well)
        if self.prepare_tables:
            for mc_set in self.nue.keys():
                self.nue[mc_set]['CC']['ebin'] = np.digitize(np.log10(self.nue[mc_set]['CC']['energy']),
                                                             self.oscCalc.prob_maps['eedges'])-1
                self.nue[mc_set]['CC']['zbin'] = np.digitize(np.cos(self.nue[mc_set]['CC']['zenith']),
                                                             self.oscCalc.prob_maps['zedges'])-1
                self.numu[mc_set]['CC']['ebin'] = np.digitize(np.log10(self.numu[mc_set]['CC']['energy']),
                                                             self.oscCalc.prob_maps['eedges'])-1
                self.numu[mc_set]['CC']['zbin'] = np.digitize(np.cos(self.numu[mc_set]['CC']['zenith']),
                                                             self.oscCalc.prob_maps['zedges'])-1
                self.nutau[mc_set]['CC']['ebin'] = np.digitize(np.log10(self.nutau[mc_set]['CC']['energy']),
                                                             self.oscCalc.prob_maps['eedges'])-1
                self.nutau[mc_set]['CC']['zbin'] = np.digitize(np.cos(self.nutau[mc_set]['CC']['zenith']),
                                                             self.oscCalc.prob_maps['zedges'])-1

        # NOTE - I don't remember if this is still necessary.
        self.oscCalc = None # Restarting oscCalc
        print 'dataLoader: MC loaded successfully!'

        return

    def loadAtmMuData(self, use_kde=False):
        '''
        Loads the templates for atmospheric muon background. 
        Only works if reconstructed values are given for the energy and zenith observables.
        '''

        self.atmmu_histo = {}

        # Atm. muons from DATA
        atmmu_data = {}
        for one_obs in self.observables:
            atmmu_data[one_obs] = np.zeros(1)

        for one_file in self.user.atmmu_data_files:
            atmmu_file = self.addCuts(pickle.load(open(one_file)), use_kde)
            for one_obs in self.observables:
                atmmu_data[one_obs] = np.concatenate((atmmu_data[one_obs], atmmu_file[one_obs]))
        self.atmmu_template = {}
        self.user.atmmu_data_files.sort()
        for one_file in self.user.atmmu_data_files:
            atmmu_file = self.addCuts(pickle.load(open(one_file)), use_kde)
            for key in atmmu_file.keys():
                if not key in self.atmmu_template.keys():
                    self.atmmu_template[key]=deepcopy(atmmu_file[key])
                else:
                    self.atmmu_template[key]= deepcopy(np.concatenate((self.atmmu_template[key], atmmu_file[key])))

        if not use_kde or len(atmmu_data[(atmmu_data.keys())[0]]) < 2:             
            self.atmmu_histo['data'], edges =  np.histogramdd([atmmu_data[x] for x in self.hist_obs], 
                                                          self.bin_edges)

        else: 
            saveFile = os.path.join(self.local_settings.kde_database,'atmmuData_')
            self.atmmu_histo['data'], edges = self.smoothByKDE( [atmmu_data[x] for x in self.hist_obs],
                                                                self.bin_edges,  
                                                                saveFile=saveFile, 
                                                                #matchKeys=matchKeys, 
                                                                #truth=truth,
                                                                mirroring=True, verbose=self.verbose )

        
        print 'dataLoader: Muon background templates loaded successfully!'
        return
    

    def loadAtmMuDataAux(self, use_kde=False):
        '''
        Loads the templates for atmospheric muon background for the MSU style of doing stuff
        Only works if reconstructed values are given for the energy and zenith observables.
        '''

        # Atm. muons from DATA
        atmmu_data = {}
        for one_obs in self.observables:
            atmmu_data[one_obs] = np.zeros(1)

        for one_file in self.user.atmmu_data_files_aux:
            atmmu_file = self.addCuts(pickle.load(open(one_file)), use_kde)
            for one_obs in self.observables:
                atmmu_data[one_obs] = np.concatenate((atmmu_data[one_obs], atmmu_file[one_obs]))
        self.atmmu_template_aux = {}
        self.user.atmmu_data_files_aux.sort()
        for one_file in self.user.atmmu_data_files_aux:
            atmmu_file = self.addCuts(pickle.load(open(one_file)), use_kde)
            for key in atmmu_file.keys():
                if not key in self.atmmu_template_aux.keys():
                    self.atmmu_template_aux[key]=deepcopy(atmmu_file[key])
                else:
                    self.atmmu_template_aux[key]= deepcopy(np.concatenate((self.atmmu_template_aux[key], 
                                                                           atmmu_file[key])))

        if not use_kde or len(atmmu_data[(atmmu_data.keys())[0]]) < 2:             
            self.atmmu_histo['data_aux'], edges =  np.histogramdd([atmmu_data[x] for x in self.hist_obs], 
                                                          self.bin_edges)

        else: 
            saveFile = os.path.join(self.local_settings.kde_database,'atmmuData_aux_')
            self.atmmu_histo['data_aux'], edges = self.smoothByKDE( [atmmu_data[x] for x in self.hist_obs],
                                                                self.bin_edges,  
                                                                saveFile=saveFile, 
                                                                #matchKeys=matchKeys, 
                                                                #truth=truth,
                                                                mirroring=True, verbose=self.verbose )

        # Rescale the aux template so it has the same integral as the one used for fitting
        if np.sum(self.atmmu_histo['data_aux']) != 0:
            self.atmmu_histo['data_aux'] *= np.sum(self.atmmu_histo['data'])/np.sum(self.atmmu_histo['data_aux'])

        
        print 'dataLoader: Muon background templates loaded successfully!'
        return


    def loadNoiseBackground(self, use_kde=False, truth={}, quantile_mc = (0., 1.)):
        '''
        Loads the templates for pure-noise simulation
        Only works if reconstructed values are given for the energy and zenith.
        '''
        # Pure noise MC
        pure_noise = {"weight":np.zeros(1)} 
        for one_obs in self.observables:
            pure_noise[one_obs] = np.zeros(1)
        for one_file in self.user.pure_noise_files:
            pure_noise_file = self.addCuts(pickle.load(open(one_file)), use_kde, quantile_mc)
            for one_obs in self.observables:
                pure_noise[one_obs] = np.concatenate((pure_noise[one_obs], pure_noise_file[one_obs]))
            pure_noise["weight"] = np.concatenate((pure_noise["weight"], pure_noise_file["weight"]))

        if not use_kde or len(pure_noise[(pure_noise.keys())[0]]) < 2: 
            self.pure_noise, edges =  np.histogramdd([pure_noise[x] for x in self.hist_obs],
                                                     self.bin_edges, weights=pure_noise['weight'])
            self.unweighted_histograms['noise_histo'], edges = \
                                                    np.histogramdd([pure_noise[x] for x in self.hist_obs], 
                                                                   self.bin_edges)
        else: 
            saveFile = os.path.join(self.local_settings.kde_database,'pureNoise_')
            self.pure_noise, edges = self.smoothByKDE(  [pure_noise[x] for x in self.hist_obs], self.bin_edges, 
                                                        saveFile=saveFile, 
                                                        # truth=truth, 
                                                        mirroring=True, verbose=self.verbose )
            # Unweighted histograms do not get any KDE business
            self.unweighted_histograms['noise_histo'], edges = \
                                                np.histogramdd([pure_noise[x] for x in self.hist_obs],
                                                self.bin_edges)
    def loadData(self, year):
        '''
        Load the experimental data. Specify the year as a string in a dictionary:
        loadData( {'year':'2011'})
        '''
        print 'dataLoader: Loading experimental data! Year ', year
        self.alldata= {}

        if year == '2011':
            data2011 = self.addCuts(pickle.load(open(self.user.data['2011'])))
            data = data2011
            self.alldata["2011"] = deepcopy(data2011)
        elif year == '2012':
            data2012 = self.addCuts(pickle.load(open(self.user.data['2012'])))
            data = data2012
            self.alldata["2012"] = deepcopy(data2012)
        elif year == '2013':
            data2013 =self.addCuts( pickle.load(open(self.user.data['2013'])))
            data = data2013
            self.alldata["2013"] = deepcopy(data2013)
        elif year == 'all':
            data2011 = self.addCuts(pickle.load(open(self.user.data['2011'])))
            data2012 = self.addCuts(pickle.load(open(self.user.data['2012'])))
            data2013 = self.addCuts(pickle.load(open(self.user.data['2013'])))
            data = {}
            self.alldata["2011"] = deepcopy(data2011)
            self.alldata["2012"] = deepcopy(data2012)
            self.alldata["2013"] = deepcopy(data2013)
            for one_obs in self.observables:
                data[one_obs] = np.concatenate((data2011[one_obs], data2012[one_obs], data2013[one_obs]))

        exp_histo, edges = np.histogramdd([data[x] for x in self.hist_obs], self.bin_edges)

        print 'dataLoader: Data loaded successfully!', np.sum(exp_histo) , ' events in final histogram.'
        return exp_histo




    def loadOriginalMCasData(self, settings, statistical_fluctuations = False):
        '''
        Creates the typical year of data as given by the MC weights.
        '''

        dsettings = deepcopy(self.default_data_settings)
        dsettings.update(settings)

        print 'dataLoader: Loading original MC as data'
        all_nu = {}
        for one_obs in self.observables:
            all_nu[one_obs] = np.concatenate((self.nue['baseline']['CC'][one_obs],  
                                              self.nue['baseline']['NC'][one_obs],  
                                              self.numu['baseline']['CC'][one_obs],   
                                              self.numu['baseline']['NC'][one_obs],   
                                              self.nutau['baseline']['CC'][one_obs], 
                                              ))
        neutrino_histo, edges = np.histogramdd([all_nu[x] for x in self.hist_obs],
                                              self.bin_edges,
                                              weights = np.concatenate((self.nue['baseline']['CC']['weight'],  
                                                                        self.nue['baseline']['NC']['weight'],
                                                                        self.numu['baseline']['CC']['weight'],
                                                                        self.numu['baseline']['NC']['weight'],
                                                                        self.nutau['baseline']['CC']['weight'], 
                                                                        )))
        exp_histo = neutrino_histo*dsettings['norm_nu']*self.sec2years
        total_neutrinos = np.sum(exp_histo)

        use_atmmu_fraction = False
        use_noise_fraction = False


	if dsettings['norm_atmmu'] == None and dsettings['atmmu_f'] == None:
	    print ' You need to set either atmmu_f or norm_atmmu. Go back to your data settings and do this.'
            sys.exit()
	if dsettings['norm_noise'] == None and dsettings['noise_f'] == None:
	    print ' You need to set either noise_f or norm_noise. Go back to your data settings and do this.'
            sys.exit()

        # Only use fractions if a normalization was not given
        if dsettings['norm_atmmu'] == None and dsettings['atmmu_f'] != False:
            print 'dataLoader: using the atmosperhic muon fraction given (not normalizations)'
            use_atmmu_fraction = True
        if dsettings['norm_noise'] == None and dsettings['noise_f'] != False:
            print 'dataLoader: using the pure noise fraction given (not normalizations)'
            use_noise_fraction = True

        # If either fraction is used, you need to know the total number of neutrinos/muons/noise produced
        if use_noise_fraction or use_atmmu_fraction:
            nu, mu, noise = self.getSingleHistos(dsettings)
            # Convert to events per year
            if dsettings['atmmu_template'] != 'data':
                mu *= self.sec2years

        # If fractions are used for everything, calculate norms
        if use_atmmu_fraction and use_noise_fraction:
            total_events = total_neutrinos*1./(1-dsettings['atmmu_f']-dsettings['noise_f'])
            dsettings['norm_atmmu'] = total_events*dsettings['atmmu_f']/np.sum(mu)
            dsettings['norm_noise'] = total_events*dsettings['noise_f']/np.sum(self.pure_noise)

        # If fractions are used only for atmmu, calculate their norm
        if use_atmmu_fraction and not use_noise_fraction:
            total_events = (total_neutrinos+dsettings['norm_noise']*np.sum(self.pure_noise))/\
                (1-dsettings['atmmu_f'])
            dsettings['norm_atmmu'] = total_events*dsettings['atmmu_f']/np.sum(mu)

        # If fractions are used only for noise, calculate their norm
        if not use_atmmu_fraction and  use_noise_fraction:
            total_events = (total_neutrinos + dsettings['norm_atmmu']*np.sum(mu))/\
                (1-dsettings['noise_f'])
            dsettings['norm_noise'] = total_events*dsettings['noise_f']/np.sum(self.pure_noise)

        print 'dataLoader: atmmu_norm ', dsettings['norm_atmmu'], ', noise_norm', dsettings['norm_noise'] 

        exp_histo += self.atmmu_histo[dsettings['atmmu_template']]*dsettings['norm_atmmu']
        exp_histo += self.pure_noise*dsettings['norm_noise']
        if statistical_fluctuations:
            exp_histo = addStatisticalFluctuations(exp_histo)

        if self.verbose:
            print 'dataLoader: Data as MC - unweighted neutrino events ', self.countMCevents()
            print 'dataLoader: Data as MC - weighted events in sample ', np.sum(exp_histo)

        return exp_histo
    #
    def normToFrac(self,settings, norm_arm, norm_noise): 
            dsettings       = deepcopy(self.default_data_settings) 
            dsettings.update(settings) 
            
            nu, mu, noise   = self.getSingleHistos(dsettings) 
            
            total_neutrinos = np.sum(nu)*dsettings['norm_nu']*self.sec2years 
            total_mu        = dsettings["norm_atmmu"] * np.sum(mu)
            if dsettings['atmmu_template'] != 'data':
                total_mu *= self.sec2years
            total_noise     = dsettings["norm_noise"] * np.sum(self.pure_noise) 
            
            print "Total: ", total_neutrinos, " ", total_mu, " ", total_noise 
            
            atmmu_f         = total_mu *1.0 / (total_neutrinos + total_mu + total_noise)  
            noise_f         = total_noise *1.0 / (total_neutrinos + total_mu + total_noise)  
            
            print "Converted norms to fractions: ", atmmu_f, noise_f 
            return atmmu_f, noise_f 
    
    def loadMCasData(self, settings, 
                     statistical_fluctuations = False):
        '''
        Creates a typical year of data using the settings received. Statistical fluctuations can be added.
        '''

        dsettings = deepcopy(self.default_data_settings)
        dsettings.update(settings)
        if self.verbose:
            print 'dataLoader: Loading MC as data. Calculating weights with the parameters that follow:'
            print dsettings

        use_atmmu_fraction = False
        use_noise_fraction = False

	if dsettings['norm_atmmu'] == None and dsettings['atmmu_f'] == None:
	    print ' You need to set either atmmu_f or norm_atmmu. Go back to your data settings and do this.'
            sys.exit()
	if dsettings['norm_noise'] == None and dsettings['noise_f'] == None:
	    print ' You need to set either noise_f or norm_noise. Go back to your data settings and do this.'
            sys.exit()

        # Only use fractions if a normalization was not given
        if dsettings['norm_atmmu'] == None and dsettings['atmmu_f'] != False:
            print 'dataLoader: using the atmosperhic muon fraction given (not normalizations)'
            use_atmmu_fraction = True
        if dsettings['norm_noise'] == None and dsettings['noise_f'] != False:
            print 'dataLoader: using the pure noise fraction given (not normalizations)'
            use_noise_fraction = True

        # If either fraction is used, you need to know the total number of neutrinos expected
        if use_noise_fraction or use_atmmu_fraction:
            nu, mu, noise = self.getSingleHistos(dsettings)
            total_neutrinos = np.sum(nu)*dsettings['norm_nu']*self.sec2years
            # Convert to events per year
            if dsettings['atmmu_template'] != 'data':
                mu *= self.sec2years

        # If fractions are used for everything, calculate norms
        if use_atmmu_fraction and use_noise_fraction:
            total_events = total_neutrinos*1./(1-dsettings['atmmu_f']-dsettings['noise_f'])
            dsettings['norm_atmmu'] = total_events*dsettings['atmmu_f']/np.sum(mu)
            dsettings['norm_noise'] = total_events*dsettings['noise_f']/np.sum(self.pure_noise)

        # If fractions are used only for atmmu, calculate their norm
        if use_atmmu_fraction and not use_noise_fraction:
            total_events = (total_neutrinos+dsettings['norm_noise']*np.sum(self.pure_noise))/\
                (1-dsettings['atmmu_f'])
            dsettings['norm_atmmu'] = total_events*dsettings['atmmu_f']/np.sum(mu)

        # If fractions are used only for noise, calculate their norm
        if not use_atmmu_fraction and  use_noise_fraction:
            total_events = (total_neutrinos + dsettings['norm_atmmu']*np.sum(mu))/\
                (1-dsettings['noise_f'])
            dsettings['norm_noise'] = total_events*dsettings['noise_f']/np.sum(self.pure_noise)

        print 'dataLoader: atmmu_norm ', dsettings['norm_atmmu'], ', noise_norm', dsettings['norm_noise'] 

        exp_histo = self.getRefHisto(params = dsettings)


        if statistical_fluctuations:
            exp_histo = addStatisticalFluctuations(exp_histo)

        if self.verbose:
            print 'dataLoader: Data as MC - unweighted neutrino events ', self.countMCevents()
            print 'dataLoader: Data as MC - weighted events in sample ', np.sum(exp_histo)

        return exp_histo

        
    def getFullNuHistogram(self, inparams = {}, weight_power = 1., use_kde=False, return_weights = False):
        '''
        Obtains the weighted MC histograms, and applies detector systematics on the histograms.
        This is the function that is finally called to obtain the MC expectation.
        '''

        params = deepcopy(self.default_data_settings)
        params.update(inparams)

        # This field needs to be declared. The value is irrelevant.
        params['baseline_correction'] = None

        # Obtaining the baseline MC
        nu_histos = self.getNuHistogram(params, weight_power=weight_power, 
                                        use_kde=use_kde, return_weights = return_weights)

        # Applying the detector-related systematics
        if params['add_detector_systematics']:

            # JP TODO - verify detailed_detsys still work (without KDEs, otherwise it's too slow)
            if self.detailed_detsys:
                redo_systematics = False
                if len(self.sysvar_params) == 0:
                    redo_systematics = True
                else:
                    for one_key in self.sysvar_keys:
                        if params[one_key] != self.sysvar_params[one_key]:
                            redo_systematics = True
                            break
                if redo_systematics:      
                    print 'Redoing systematics'
                    self.loadSystematicVariations(params = params, use_kde=use_kde)

            # Apply the bin-wise correction in a flavor-wise way
            # Baseline corrections are applied as another systematic
            for particle in self.nu_keys:
                for i in range(nu_histos[particle+'_histo'].shape[0]):
                    for j in range(nu_histos[particle+'_histo'].shape[1]):
                        for k in range(nu_histos[particle+'_histo'].shape[2]):
                            bin_factor = 1.
                            for systKey in self.sysfcn[particle].keys():
                                bin_factor *= self.sysfcn[particle][systKey][i,j,k](params[systKey])
                            nu_histos[particle+'_histo'][i,j,k] *= bin_factor**weight_power
                            self.last_detsys_factors[particle][i+1,j+1,k+1] = bin_factor**weight_power

        return nu_histos
        
                        
    def getNuHistogram(self, inparams = {}, weight_power = 1., use_kde=False, return_weights = False):
        '''
        Obtains the weighted MC histograms WITHOUT detector systematics!
        '''

        params = deepcopy(self.default_data_settings)
        params.update(inparams)
                       
        parameters_changed = False
        # If the parameters change, get the calculator again
        if (self.oscCalc == None or
            self.oscCalc.oscMode != params['oscMode'] or 
            self.oscCalc.tables != params['oscTables']):
            self.oscCalc = oscCalc.OscCalc(oscMode  = params['oscMode'],
                                           doTables = params['oscTables'],
                                           nbins = self.iniDict['table_nbins'])
            parameters_changed = True


        if (self.oscCalc.dm31 != params['dm31'] or 
            self.oscCalc.theta23 != params['theta23'] or 
            self.oscCalc.theta13 != params['theta13'] or
            self.oscCalc.mix_angle != params['mix_angle'] or 
            self.oscCalc.dm41 != params['dm41'] or
            self.oscCalc.theta24 != params['theta24'] or 
            self.oscCalc.theta34 != params['theta34'] ):

            self.oscCalc.setParameters(dm31=params['dm31'],
                                       theta23=params['theta23'], 
                                       theta13=params['theta13'],
                                       mix_angle=params['mix_angle'],
                                       dm41=params['dm41'], 
                                       theta24=params['theta24'],
                                       theta34=params['theta34']) 
            parameters_changed = True



        nue_cc = self.nue[params['simulation']]['CC']
        numu_cc = self.numu[params['simulation']]['CC']
        nutau_cc = self.nutau[params['simulation']]['CC']

        #print 'Test', len(nue_cc['weight_mu']), len(numu_cc['weight_e']), len(nutau_cc['weight_mu'])
        # Calculation of oscillation probabilities
        # If the parameters did not change, and the simulation is the same, do nothing
        # If the parameters changed, then re-calculate
        # If the simulation changed, then re-calculate
        # Define selection of DIS, QE, RES here
        if self.oscMC != params['simulation'] or parameters_changed:

            self.oscMC = params['simulation']

            if params['oscTables']:
                energy_key = 'ebin'
                zenith_key = 'zbin'
            else:
                energy_key = 'energy'
                zenith_key = 'zenith'                
            #print nue_cc.keys(), nue_cc['ptype']
            self.P_nuX_to_nue   = self.oscCalc.GetProb(nue_cc[energy_key], 
                                                       nue_cc[zenith_key], 
                                                       nue_cc['ptype'])
            self.P_nuX_to_numu  = self.oscCalc.GetProb(numu_cc[energy_key], 
                                                       numu_cc[zenith_key], 
                                                       numu_cc['ptype'])
            self.P_nuX_to_nutau = self.oscCalc.GetProb(nutau_cc[energy_key], 
                                                       nutau_cc[zenith_key], 
                                                       nutau_cc['ptype'])

            self.nuscatt = {'nue':  {'DIS' : self.nue[self.oscMC]['CC']['scattering'] == 1,
                                     'NCDIS':self.nue[self.oscMC]['NC']['scattering'] == 1,
                                     'RES' : self.nue[self.oscMC]['CC']['scattering'] == 2,
                                     'QEL' : self.nue[self.oscMC]['CC']['scattering'] == 3},
                            'numu': {'DIS' : self.numu[self.oscMC]['CC']['scattering'] == 1,
                                     'NCDIS':self.numu[self.oscMC]['NC']['scattering'] == 1,
                                     'RES' : self.numu[self.oscMC]['CC']['scattering'] == 2,
                                     'QEL' : self.numu[self.oscMC]['CC']['scattering'] == 3},
                            'nutau':{'DIS' : self.nutau[self.oscMC]['CC']['scattering'] == 1,
                                     'NCDIS':self.nutau[self.oscMC]['NC']['scattering'] == 1,
                                     'RES' : self.nutau[self.oscMC]['CC']['scattering'] == 2,
                                     'QEL' : self.nutau[self.oscMC]['CC']['scattering'] == 3},
                            }


        # Calculation of event weights, including flux and cross-sections

        # All flux calcualtions are 
        # Flux = (original_flux*spectral_index_mod*nu/nubar_ratio)*nue/numu_ratio (for nuE source only)
        # The parameters that modify the flux were derived from the Barr paper

        # NuE produced as NuE
        nue_nue_ccflux    = nue_cc[self.we]*params['norm_e']*nue_cc['energy']**params['gamma']*\
                            sf.modRatioNuEBar(nue_cc, params['nu_nubar'], params['nubar_ratio'])*\
                            sf.modRatioUpHor_NuE(nue_cc, params['uphor_ratio'])


        # NuE produced as NuMu
        nue_numu_ccflux   = (nue_cc[self.wmu]+params['nu_pi_scale']*nue_cc[self.wmu_pi])*\
                            nue_cc['energy']**params['gamma']*\
                            sf.modRatioNuMuBar(nue_cc, params['nu_nubar'], params['nubar_ratio'])*\
                            sf.modRatioUpHor_NuMu(nue_cc, params['uphor_ratio'])

        # NuMu produced as NuE
        numu_nue_ccflux   = numu_cc[self.we]*params['norm_e']*numu_cc['energy']**params['gamma']*\
                            sf.modRatioNuEBar(numu_cc, params['nu_nubar'], params['nubar_ratio'])*\
                            sf.modRatioUpHor_NuE(numu_cc, params['uphor_ratio'])
        
        # NuMu produced as NuMu
        numu_numu_ccflux  = (numu_cc[self.wmu]+params['nu_pi_scale']*numu_cc[self.wmu_pi])*\
                            numu_cc['energy']**params['gamma']*\
                            sf.modRatioNuMuBar(numu_cc, params['nu_nubar'], params['nubar_ratio'])*\
                            sf.modRatioUpHor_NuMu(numu_cc, params['uphor_ratio'])

        # NuTau produced as NuE
        nutau_nue_ccflux    = nutau_cc[self.we]*params['norm_e']*nutau_cc['energy']**params['gamma']*\
                              sf.modRatioNuEBar(nutau_cc, params['nu_nubar'], params['nubar_ratio'])*\
                              sf.modRatioUpHor_NuE(nutau_cc, params['uphor_ratio'])
        
        # NuTau produced as NuMu
        nutau_numu_ccflux   = (nutau_cc[self.wmu]+params['nu_pi_scale']*nutau_cc[self.wmu_pi])*\
                              nutau_cc['energy']**params['gamma']*\
                              sf.modRatioNuMuBar(nutau_cc, params['nu_nubar'], params['nubar_ratio'])*\
                              sf.modRatioUpHor_NuMu(nutau_cc, params['uphor_ratio'])
        

        # NC - oscillations not calculated. To calculate for NC, just do the previous steps as well.
        nue_nc = self.nue[self.oscMC]['NC']
        numu_nc = self.numu[self.oscMC]['NC']
 
        # Calculating NC oscillations to sterile state
        if self.oscMC != params['simulation'] or parameters_changed:
            if params['oscMode'].lower() == "globes_sterile":
                self.P_nue_to_nus = self.oscCalc.GetProbToNuS( nue_nc['energy'],
                                                               nue_nc['zenith'],nue_nc['ptype'])
                self.P_numu_to_nus= self.oscCalc.GetProbToNuS(numu_nc['energy'],
                                                              numu_nc['zenith'],numu_nc['ptype'])  

        nue_ncweight  = nue_nc[self.we]*params['norm_e']*nue_nc['energy']**params['gamma']*\
                        sf.modRatioNuEBar(nue_nc, params['nu_nubar'], params['nubar_ratio'])*\
                        sf.modRatioUpHor_NuE(nue_nc, params['uphor_ratio'])


        numu_ncweight  = (numu_nc[self.wmu]+params['nu_pi_scale']*numu_nc[self.wmu_pi])*\
                         numu_nc['energy']**params['gamma']*\
                         sf.modRatioNuMuBar(numu_nc, params['nu_nubar'], params['nubar_ratio'])*\
                         sf.modRatioUpHor_NuMu(numu_nc, params['uphor_ratio'])

        # Applying HONDA16 correction to fluxes which are pre-2016
        # The correction is applied to all flavors, and is only energy-dependent
        # This section could be moved to the data loading. I placed it here temporarily for testing!
        if params['correct_honda_flux']:
            nuecc_flux_corr  = sf.honda_fix(nue_cc['energy'])
            numucc_flux_corr = sf.honda_fix(numu_cc['energy'])
            nutaucc_flux_corr= sf.honda_fix(nutau_cc['energy'])
            nuenc_flux_corr  = sf.honda_fix(nue_nc['energy'])
            numunc_flux_corr = sf.honda_fix(numu_nc['energy'])

            nue_nue_ccflux    *= nuecc_flux_corr
            nue_numu_ccflux   *= nuecc_flux_corr
            numu_nue_ccflux   *= numucc_flux_corr
            numu_numu_ccflux  *= numucc_flux_corr
            nutau_nue_ccflux  *= nutaucc_flux_corr
            nutau_numu_ccflux *= nutaucc_flux_corr
            nue_ncweight      *= nuenc_flux_corr
            numu_ncweight     *= numunc_flux_corr



        # Applying oscillations to the unoscillated flux
        # NC to sterile first
        if params['oscMode'].lower() == "globes_sterile": 
            nue_ncweight = (1.0 - self.P_nue_to_nus[0])*nue_ncweight
            numu_ncweight = (1.0 - self.P_numu_to_nus[0])*numu_ncweight
        if params['dm31'] != 0:
            nue_ccweight  =self.P_nuX_to_nue[:,0]*nue_nue_ccflux+self.P_nuX_to_nue[:,1]*nue_numu_ccflux
            numu_ccweight =self.P_nuX_to_numu[:,0]*numu_nue_ccflux+self.P_nuX_to_numu[:,1]*numu_numu_ccflux
            nutau_ccweight=(self.P_nuX_to_nutau[:,0]*nutau_nue_ccflux+self.P_nuX_to_nutau[:,1]*nutau_numu_ccflux)*\
                params['norm_tau']
        else:
            nue_ccweight  =nue_nue_ccflux
            numu_ccweight =numu_numu_ccflux
            nutau_ccweight=np.zeros_like(nutau_nue_ccflux)


        # CROSS SECTION UNCERTAINTIES
        # Retrieving the changes due to the axial mass variations
        if params['ma_variations']:
            nue_ccweight[self.nuscatt['nue']['RES']]   *= \
                sf.axialMassVar(coeff = self.nue[self.oscMC]['CC']['ma_res'][self.nuscatt['nue']['RES'],:],
                                Ma = params['axm_res'])
            nue_ccweight[self.nuscatt['nue']['QEL']]   *= \
                sf.axialMassVar(coeff = self.nue[self.oscMC]['CC']['ma_qe'][self.nuscatt['nue']['QEL'],:],
                                Ma = params['axm_qe'])

            numu_ccweight[self.nuscatt['numu']['RES']]   *= \
                sf.axialMassVar(coeff = self.numu[self.oscMC]['CC']['ma_res'][self.nuscatt['numu']['RES'],:],
                                Ma = params['axm_res'])
            numu_ccweight[self.nuscatt['numu']['QEL']]   *= \
                sf.axialMassVar(coeff = self.numu[self.oscMC]['CC']['ma_qe'][self.nuscatt['numu']['QEL'],:],
                                Ma = params['axm_qe'])

            nutau_ccweight[self.nuscatt['nutau']['RES']]   *= \
                sf.axialMassVar(coeff = self.nutau[self.oscMC]['CC']['ma_res'][self.nuscatt['nutau']['RES'],:],
                                Ma = params['axm_res'])
            nutau_ccweight[self.nuscatt['nutau']['QEL']]   *= \
                sf.axialMassVar(coeff = self.nutau[self.oscMC]['CC']['ma_qe'][self.nuscatt['nutau']['QEL'],:],
                                Ma = params['axm_qe'])


        # Applying DIS cross section uncertainty ('scattering' key should be 1 in pickle files)
        if params['dis_reweight']:
            nue_ccweight[self.nuscatt['nue']['DIS']] *= \
                sf.tkDISreweight(a=params['dis_xa'], b=params['dis_xb'], 
                                 bjorken_x=self.nue[self.oscMC]['CC']['GENIE_x'][self.nuscatt['nue']['DIS']])
            nue_ncweight[self.nuscatt['nue']['NCDIS']] *= \
                sf.tkDISreweight(a=params['dis_xa'], b=params['dis_xb'], 
                                 bjorken_x=self.nue[self.oscMC]['NC']['GENIE_x'][self.nuscatt['nue']['NCDIS']])
            numu_ccweight[self.nuscatt['numu']['DIS']] *= \
                sf.tkDISreweight(a=params['dis_xa'], b=params['dis_xb'], 
                                 bjorken_x=self.numu[self.oscMC]['CC']['GENIE_x'][self.nuscatt['numu']['DIS']])
            numu_ncweight[self.nuscatt['numu']['NCDIS']] *= \
                sf.tkDISreweight(a=params['dis_xa'], b=params['dis_xb'], 
                                 bjorken_x=self.numu[self.oscMC]['NC']['GENIE_x'][self.nuscatt['numu']['NCDIS']])

            nutau_ccweight[self.nuscatt['nutau']['DIS']] *= \
                sf.tkDISreweight(a=params['dis_xa'], b=params['dis_xb'], 
                                 bjorken_x=self.nutau[self.oscMC]['CC']['GENIE_x'][self.nuscatt['nutau']['DIS']])



        # Merging NC into one component
        nu_nc = {}
        for one_obs in self.observables:
            nu_nc[one_obs] = np.concatenate((nue_nc[one_obs], numu_nc[one_obs]))
        nu_nc_weight    = params['norm_nc']*np.concatenate((nue_ncweight, numu_ncweight))

        # Do not modify the weights after this point
        if return_weights:
            self.nue[self.oscMC]['CC']['out_weight']   = nue_ccweight
            self.nue[self.oscMC]['NC']['out_weight']   = params['norm_nc']*nue_ncweight
            self.numu[self.oscMC]['CC']['out_weight']  = numu_ccweight
            self.numu[self.oscMC]['NC']['out_weight']  = params['norm_nc']*numu_ncweight
            self.nutau[self.oscMC]['CC']['out_weight'] = nutau_ccweight


        # Apply change in PID (it is undone after the histogram is calculated)
        nue_cc[self.observables[-1]]   += params['pid_bias']
        numu_cc[self.observables[-1]]  += params['pid_bias']
        nutau_cc[self.observables[-1]] += params['pid_bias']
        nu_nc[self.observables[-1]]   += params['pid_bias']

        # Placing all the information in histograms
        if not use_kde:
            numu_histo, edges = np.histogramdd([numu_cc[x] for x in self.hist_obs],
                                           self.bin_edges, weights=numu_ccweight**weight_power)
            nue_histo, edges  = np.histogramdd([nue_cc[x] for x in self.hist_obs],
                                           self.bin_edges, weights = nue_ccweight**weight_power)
            nutau_histo,edges = np.histogramdd([nutau_cc[x] for x in self.hist_obs],
                                           self.bin_edges, weights = nutau_ccweight**weight_power)
            nc_histo, edges   = np.histogramdd([nu_nc[x] for x in self.hist_obs],
                                           self.bin_edges,
                                           weights=nu_nc_weight**weight_power)
        else: 
            if self.oscMC in ("baseline", "Baseline"):  tFile = self.baseline_key 
            else:                                       tFile = self.oscMC 

            matchKeys   = [ 'pid_bias', 'oscMode', 'dm31', 'theta23', 'theta13', 
                            'uphor_ratio', 'nu_nubar', 'nubar_ratio', 'oscTables', 
                            'simulation', 'norm_e', 'gamma', "weight_power" ] 

            tTruth            = deepcopy(  params  ) 
            tTruth.update(   { "weight_power"        : weight_power          }  ) 
            tTruth.update(   { "numu_weight_sum"     : np.sum(numu_ccweight) }  ) 
          
            saveFile = os.path.join(self.local_settings.kde_database, self.oscMC + '_' + "%i" % weight_power + '_')
            numu_histo, edges = self.smoothByKDE(   [numu_cc[x] for x in self.hist_obs], 
                                                    self.bin_edges, 
                                                    weights=numu_ccweight**weight_power, 
                                                    truth=tTruth,
                                                    saveFile = saveFile+'numu.pckl',
                                                    matchKeys=matchKeys, 
                                                    verbose=self.verbose )
            nue_histo, edges  = self.smoothByKDE(   [nue_cc[x] for x in self.hist_obs], 
                                                    self.bin_edges, 
                                                    weights = nue_ccweight**weight_power, 
                                                    truth=tTruth, 
                                                    saveFile = saveFile+'nue.pckl',
                                                    matchKeys=matchKeys, 
                                                    verbose=self.verbose )
            nutau_histo,edges = self.smoothByKDE(   [nutau_cc[x] for x in self.hist_obs], 
                                                    self.bin_edges, 
                                                    weights = nutau_ccweight**weight_power, 
                                                    truth=tTruth, 
                                                    saveFile = saveFile+'nutau.pckl',
                                                    matchKeys=matchKeys, 
                                                    verbose=self.verbose )
            nc_histo, edges   = self.smoothByKDE(   [nu_nc[x] for x in self.hist_obs], 
                                                    self.bin_edges, 
                                                    weights=nu_nc_weight**weight_power, 
                                                    truth=tTruth, 
                                                    saveFile = saveFile+'nc.pckl',
                                                    matchKeys=matchKeys, 
                                                    fileSuffix="_NC", 
                                                    verbose=self.verbose )
        # Undo the PID change
        nue_cc[self.hist_obs[-1]]   -= params['pid_bias']
        numu_cc[self.hist_obs[-1]]  -= params['pid_bias']
        nutau_cc[self.hist_obs[-1]] -= params['pid_bias']
        nu_nc[self.hist_obs[-1]]    -= params['pid_bias']


        neutrino_histograms = {'numu_histo':numu_histo, 'nue_histo':nue_histo, 
                               'nutau_histo':nutau_histo, 'nc_histo':nc_histo}           

        return neutrino_histograms


    def getAtmMuHistogram(self, inparams = {}, weight_power = 1., use_kde = False, return_weights = False):
        '''
        Obtain the weighted muon MC histograms WITHOUT detector systematics
        '''

        params = deepcopy(self.default_data_settings)
        params.update(inparams)
        
        # TODO - Right now nothing happens to the atmospheric muons. This will have to change

        # Apply the changes in the PID 
        self.atmmu[params['simulation']][self.observables[-1]] += params['pid_bias']

        # Placing all the info in a histogram
        if not use_kde:
            # JPMod  - change the weight key here. Make it dynamic
            atmmu_histo, edges = np.histogramdd([self.atmmu[params['simulation']][x] for x in self.hist_obs],
                                               self.bin_edges, 
                                               weights=self.atmmu[params['simulation']]['weight_mu']**weight_power)
        else:
            saveFile = os.path.join(self.local_settings.kde_database,
                                    params['simulation'] + '_' + "%i" % weight_power + '_')
            matchKeys = ['gamma']

            atmmu_histo, edges = self.smoothByKDE([self.atmmu[params['simulation']][x] for x in self.hist_obs],
                                        self.bin_edges, 
                                        weights=self.atmmu[params['simulation']]['weight_mu']**weight_power,
                                        truth = deepcopy(params),
                                        saveFile = saveFile + 'atmmu.pckl',
                                        matchKeys = matchKeys,
                                        verbose=self.verbose)
        

        # Undo the PID change
        self.atmmu[params['simulation']][self.observables[-1]] -= params['pid_bias']

        # JPmod - I haven't done much to the weights, so I just copy them for now. Modify this.
        if return_weights:
            self.atmmu[params['simulation']]['out_weight'] = self.atmmu[params['simulation']]['weight']

        return atmmu_histo

    def getFullAtmMuHistogram(self, inparams = {}, weight_power = 1., use_kde = False, return_weights = False):
        '''
        Obtains the weighted MC histograms, and applies detector systematics on the histograms.
        This is the function that is finally called to obtain the MC expectation.
        '''

        params = deepcopy(self.default_data_settings)
        params.update(inparams)

        # This field needs to be declared. The value is irrelevant.
        params['baseline_correction'] = None

        # Obtaining the baseline AtmMu
        atmmu_histo = self.getAtmMuHistogram(params,  weight_power = 1.,
                                             use_kde = use_kde, return_weights = return_weights)

        # Applying the detector-related systematics
        if params['add_detector_systematics']:

            if self.detailed_detsys:
                print 'dataLoader: Method not implemented yet'
            
            for i in range(atmmu_histo.shape[0]):
                for j in range(atmmu_histo.shape[1]):
                    for k in range(atmmu_histo.shape[2]):
                        bin_factor = 1.
                        for systKey in self.sysfcn['atmmu'].keys():
                            bin_factor *= self.sysfcn['atmmu'][systKey][i,j,k](params[systKey])
                        atmmu_histo[i,j,k] *= bin_factor*weight_power
                        self.last_detsys_factors['atmmu'][i+1,j+1,k+1] = bin_factor*weight_power

        return atmmu_histo


    def calculateEventParams(self, params = None , output_name = False):
        '''
        Calculate the weights event-wise for a given set of parameters
        '''

        if params == None:
            print 'dataLoader ERROR: You need to give all parameters when asking for event weights'

        if output_name:
            out_data = {}

        # These instructions will get the detsys factors and the neutrino/muon weights
        self.getFullNuHistogram(inparams = params, return_weights = True)
        self.getFullAtmMuHistogram(inparams = params, return_weights = True)

        params['norm_nu'] *= self.sec2years
        if params['atmmu_template'] != 'data':
            params['norm_atmmu'] *= self.sec2years 

      
        # Loop over the particles
        for particle, norm, pdict in zip( [ 'nue','nc','numu','nc','nutau','atmmu'],
                                          [ params['norm_nu']] *5 + [params['norm_atmmu']],
                                          [ self.nue[self.oscMC]['CC'], self.nue[self.oscMC]['NC'],
                                            self.numu[self.oscMC]['CC'], self.numu[self.oscMC]['NC'],
                                            self.nutau[self.oscMC]['CC'],self.atmmu['baseline'] ]):
           
            # Special case for no muon MC
            if not bool(pdict):
                continue

            # Find their position in each dimension
            pdict_length = pdict[self.hist_obs[0]].size
            indices = np.zeros([pdict_length, 3], dtype=int)
            
            # Digitize gives 0 or len(edges) when out of bounds. Taking advantage of it to produce zero weights.
            for i, one_obs in enumerate(self.hist_obs):
                indices[:,i] = np.digitize(pdict[one_obs], self.bin_edges[i])
            
            # Apply the correction factor - maybe there's a smarter way to evaluate an array with an array
            for iev in range(pdict_length):
                pdict['out_weight'][iev] *= self.last_detsys_factors[particle][zip(indices[iev,:])][0]
            
            pdict['out_weight'] *= norm

            if output_name:
                if out_data.has_key(particle):
                    for one_key in out_data[particle].keys():
                        out_data[particle][one_key] = np.concatenate((out_data[particle][one_key],
                                                                      pdict[one_key]))
                else:
                    out_data[particle] = deepcopy(pdict)

        if output_name:
            pickle.dump(out_data, open(os.path.join(self.user.data_dir, output_name + '.pckl'), 'w'))
        
            


    def getRefHisto(self, params = None):
        '''
        Returns a prediction for a given set of parameters, including neutrinos, muons and noise.
        The histogram is summed according to the normalizations provided.
        '''

        if params == None:
            print 'dataLoader ERROR: You need to define parameters for the reference histogram'

        neutrino_histo, atmmu_histo, noise_histo = self.getSingleHistos(params = params)

        ref_histo = (params['norm_nu']*self.sec2years*neutrino_histo + 
                     params['norm_noise']*self.pure_noise)

        # Muon MC needs to be scaled to years as well
        if params['atmmu_template'] != 'data':
            ref_histo += atmmu_histo*params['norm_atmmu']*self.sec2years
        else:
            ref_histo += atmmu_histo*params['norm_atmmu']

        self.norm_atmmu_last_call = deepcopy(params['norm_atmmu'])
        if self.verbose:
            print 'Using the following parameters for the reference histogram'
            print params
            print 'Events in ref_histo', np.sum(ref_histo)
        return ref_histo


    def getSingleHistos(self, params = None, detailed_nu = False, use_kde_bg=False, weight_power = 1.):
        '''
        Returns the neutrinos (weighted) and muons separately. The user (or fitter) must handle normalizations.
        '''
        nuh             = self.getFullNuHistogram(inparams = params, weight_power = weight_power)
        neutrino_histo  = (nuh['numu_histo'] + nuh['nue_histo'] + nuh['nc_histo'] + nuh['nutau_histo']) 

        if params['atmmu_template'] == 'mc' or params['atmmu_template'] == 'corsika':
            atmmu_histo = self.getFullAtmMuHistogram(inparams = params, use_kde=use_kde_bg, 
                                                     weight_power=weight_power)
        elif params['atmmu_template'] == 'data':
            atmmu_histo     = self.atmmu_histo[params['atmmu_template']]

        noise           = self.pure_noise
        
        if detailed_nu:
            return neutrino_histo, atmmu_histo, noise, nuh
        else:
            return neutrino_histo, atmmu_histo, noise


    def stageSystematicVariations(self):
        '''
        Defines necessary variables for parameterizing detector uncertainties.
        '''

        # Initializing the detector systematics
        self.detsys = {}
        tolerance = 1E-9

        # Verify that all the systematics requested are part of the user definition
        use_syst = []
        for systKey in self.detsys_nuspecs.keys():
            if not systKey in self.user.mc_sets.keys():
                print 'dataLoader: Requesting', systKey, 'nu systematic which is not defined in user. Removing.'
                self.detsys_nuspecs.pop(systKey, None)

        for systKey in self.detsys_muspecs.keys():
            if not systKey in self.user.atmmu_sets.keys():
                print 'dataLoader: Requesting', systKey, 'mu systematic which is not defined in user. Removing.'
                self.detsys_muspecs.pop(systKey, None)


        # Neutrinos - Some info is duplicated because it simplifies further steps
        for particle in self.nu_keys:
            self.detsys[particle] = {}
            # Get the evaluation points and the nominal MC
            for i, systKey in enumerate(self.detsys_nuspecs.keys()):
                xpoints = np.array(self.user.mc_sets[systKey].keys())
                xpoints.sort()
                nominal_diff = np.abs(xpoints - self.detsys_nuspecs[systKey][0])
                xnominal = nominal_diff.argmin()
                if nominal_diff[xnominal] > (self.detsys_nuspecs[systKey][0]*tolerance):
                    print 'dataLoader: ERROR! nominal value for ', systKey, ' nu MC set NOT loaded. Check value!'
                    sys.exit()    
                self.detsys[particle][systKey] = {'xpoints': xpoints,
                                                  'xnominal': xnominal,
                                                  'poldeg': self.detsys_nuspecs[systKey][1]
                                                  }
        # Muons - Some info is duplicated because it simplifies further steps
        particle = 'atmmu'
        self.detsys[particle] = {}
        for i, systKey in enumerate(self.detsys_muspecs.keys()):
            xpoints = np.array(self.user.atmmu_sets[systKey].keys())
            xpoints.sort()
            nominal_diff = np.abs(xpoints - self.detsys_muspecs[systKey][0])
            xnominal = nominal_diff.argmin()
            if nominal_diff[xnominal] > (self.detsys_muspecs[systKey][0]*tolerance):
                print 'dataLoader: ERROR! nominal value for ', systKey, ' mu MC set NOT loaded. Check value!'
                sys.exit()    
            self.detsys[particle][systKey] = {'xpoints': xpoints,
                                              'xnominal': xnominal,
                                              'poldeg': self.detsys_muspecs[systKey][1]
                                              }


    def loadSystematicVariations(self, params = None, use_kde=False):

        '''
        Loads/produces the histograms for the fit of systematic uncertainties. Fit done afterwards.
        Done for each neutrino interaction and for atmospheric muons.
        '''

        # Forcing the baseline simulation
        params['simulation'] = 'baseline'
        self.sysvar_params = deepcopy(params)
        
        # NEUTRINOS
        baseline = self.getNuHistogram(inparams = params, use_kde=use_kde) 
        # JP check - hardcoded errors NOT to use KDE function
        baseline_w2  = self.getNuHistogram(inparams = params, weight_power = 2., use_kde=False) 
        
        # Loop over the different sources of error
        for systKey in self.detsys[self.nu_keys[0]].keys():
            
            # Initialize empty systematic histograms
            for nu in self.nu_keys:
                self.detsys[nu][systKey]['hist_shape'] = baseline[nu+'_histo'].shape
                self.detsys[nu][systKey]['y'] = np.zeros( [len(self.detsys[nu][systKey]['xpoints'])] +
                                                          list(baseline[nu + '_histo'].shape))
                self.detsys[nu][systKey]['yw2'] = np.zeros_like(self.detsys[nu][systKey]['y'])

            # Populate the histograms
            for index, systValue in enumerate(self.detsys[self.nu_keys[0]][systKey]['xpoints']):
                if index == self.detsys[self.nu_keys[0]][systKey]['xnominal']:
                    tmp_hist_dict   = baseline
                    tmp_hist_dictw2 = baseline_w2
                else:
                    params['simulation'] = self.user.mc_sets[systKey][systValue]
                    tmp_hist_dict    = self.getNuHistogram(inparams=params, use_kde=use_kde)
                    # JP check - hardcoded errors NOT to use KDE function
                    tmp_hist_dictw2  = self.getNuHistogram(inparams=params, weight_power=2., use_kde=False)

                # Assign the histograms to the appropriate neutrino flavor
                for nu in self.nu_keys:
                    self.detsys[nu][systKey]['y'][index]   = tmp_hist_dict[nu+'_histo']
                    self.detsys[nu][systKey]['yw2'][index] = tmp_hist_dictw2[nu+'_histo']

        # If there aren't any atmospheric muons, go out here
        if len(self.user.atmmu_sets.keys()) == 0:
            print 'dataLoader: No MC for atmospheric muons'
            return

        # ATMOSPHERIC MUONS
        params['simulation'] = 'baseline'
        baseline = self.getAtmMuHistogram(inparams = params, use_kde = use_kde)
        baseline_w2 = self.getAtmMuHistogram(inparams = params, weight_power = 2., use_kde = False)

        # Loop over the different sources of error
        for systKey in self.detsys['atmmu'].keys():
            self.detsys['atmmu'][systKey]['hist_shape'] = baseline.shape
            self.detsys['atmmu'][systKey]['y'] = np.zeros( [len(self.detsys['atmmu'][systKey]['xpoints'])] +
                                                           list(baseline.shape))
            self.detsys['atmmu'][systKey]['yw2'] = np.zeros_like(self.detsys['atmmu'][systKey]['y'])

            # Populate the histograms (assign them already)
            for index, systValue in enumerate(self.detsys['atmmu'][systKey]['xpoints']):
                if index == self.detsys['atmmu'][systKey]['xnominal']:
                    self.detsys['atmmu'][systKey]['y'][index] = baseline
                    self.detsys['atmmu'][systKey]['yw2'][index] = baseline_w2
                else:
                    params['simulation'] = self.user.atmmu_sets[systKey][systValue]
                    self.detsys['atmmu'][systKey]['y'][index] = self.getAtmMuHistogram(inparams=params,
                                                                                             use_kde = use_kde)
                    self.detsys['atmmu'][systKey]['yw2'][index] = self.getAtmMuHistogram(inparams=params,
                                                                                               weight_power = 2.,
                                                                                               use_kde = False)
        return

    def buildSystematicFunctions(self):
        '''
        Receive the systematic histograms and build the systematic functions.
        Creates coefficients as an intermediate step.
        '''

        self.sysfcn = {}

        # Need to build each systematic function
        for particle in self.nu_keys+['atmmu']:
            self.sysfcn[particle] = {}
            for systKey in self.detsys[particle].keys():
                
                # Do the fit, store the coefficients
                print 'dataLoader: Systematic functions for', particle, ' - ', systKey
                self.detsys[particle][systKey]['coeffs'], self.detsys[particle][systKey]['covmatrix'] = \
                                sf.getSystematicFunctions(self.detsys[particle][systKey]['y'],
                                                          self.detsys[particle][systKey]['yw2'],
                                                          self.detsys[particle][systKey]['xpoints'],
                                                          self.detsys[particle][systKey]['xnominal'],
                                                          self.detsys[particle][systKey]['poldeg'])

                # Build the functions out of the coefficients
                self.sysfcn[particle][systKey] = np.zeros(self.detsys[particle][systKey]['hist_shape'],
                                                          dtype=object)
                # Triple loop to go over each bin
                for i in range(self.sysfcn[particle][systKey].shape[0]):
                    for j in range(self.sysfcn[particle][systKey].shape[1]):
                        for k in range(self.sysfcn[particle][systKey].shape[2]):
                            # Depending on the polynomial degree
                            if self.detsys[particle][systKey]['poldeg'] > 0:
                                self.sysfcn[particle][systKey][i,j,k] = \
                                            np.poly1d(self.detsys[particle][systKey]['coeffs'][i,j,k])
                            else:
                                self.sysfcn[particle][systKey][i,j,k] = \
                                            sf.exp_func(*self.detsys[particle][systKey]['coeffs'][i,j,k], 
                                                       x0 = self.detsys_muspecs[systKey][0])
        return

    def combineSystematics(self):

        # Initialize the factors storage (bins+2 size)
        self.last_detsys_factors = {}

        # Loop over each systematic (per particle)
        for particle in self.nu_keys+['atmmu']:

            self.last_detsys_factors[particle] = np.zeros([x.size+1 for x in self.bin_edges])
            print 'dataLoader: Combining systematics for ', particle

            if len(self.sysfcn[particle].keys()) == 0:
                print 'dataLoader: Not building systematic functions for', particle
                continue

            # Start two arrays to store the alignment
            one_syst = self.sysfcn[particle].keys()[0]
            align_sys_wsum     = np.zeros(self.detsys[particle][one_syst]['hist_shape'])
            align_sys_wcorrsum = np.zeros(self.detsys[particle][one_syst]['hist_shape'])

            for systKey in self.sysfcn[particle].keys():

                # Add the error and correction dictionaries
                self.detsys[particle][systKey]['corr'] = np.ones(self.sysfcn[particle][systKey].shape)
                self.detsys[particle][systKey]['corr_err'] = np.ones(self.sysfcn[particle][systKey].shape)


                for i in range(self.sysfcn[particle][systKey].shape[0]):
                    for j in range(self.sysfcn[particle][systKey].shape[1]):
                        for k in range(self.sysfcn[particle][systKey].shape[2]):                
                            
                            xnominal = \
                            self.detsys[particle][systKey]['xpoints'][self.detsys[particle][systKey]['xnominal']]
                            
                            # JPmod ... shamelessly copied from Martin. Don't know what he's doing
                            if self.detsys[particle][systKey]['poldeg'] > 0:
                                vec=[xnominal**s for s in range(len(self.sysfcn[particle][systKey][i,j,k].coeffs))]
                            else:
                                vec=self.sysfcn[particle][systKey][i,j,k].derive(xnominal)

                            # Calculating the error on y from fit at x = nominal
                            tProd = np.dot(np.dot(vec, self.detsys[particle][systKey]['covmatrix'][i,j,k]), vec)
                            err   = np.sqrt(np.abs(tProd))
                            # JPmod - Martin put a warning here on weird behavior of cov. matrix
                            
                            # Calculate correction for sys spline due to nominal bias
                            tCorr = 1.0/self.sysfcn[particle][systKey][i,j,k](xnominal)
                            self.detsys[particle][systKey]['corr'][i,j,k] = tCorr
                            self.detsys[particle][systKey]['corr_err'][i,j,k] = ((1.0/tCorr)**2)*err
                            
                            if self.detsys[particle][systKey]['poldeg'] > 0:
                                self.sysfcn[particle][systKey][i,j,k] *= tCorr
                            else:
                                self.sysfcn[particle][systKey][i,j,k].a *= tCorr
                                self.sysfcn[particle][systKey][i,j,k].c *= tCorr

                            # JPmod -  Didn't include the warnings of Martin so far
                            
                            # Pass information about nominal bias
                            align_sys_wsum[i,j,k]    +=   1./self.detsys[particle][systKey]['corr_err'][i,j,k]**2
                            align_sys_wcorrsum[i,j,k]+=tCorr/self.detsys[particle][systKey]['corr_err'][i,j,k]**2

            # Implementing a dummy function that always returns a fixed value - generalizes evaluations later
            # Using the last systKey defined
            self.sysfcn[particle]['baseline_correction'] = np.zeros(self.detsys[particle][systKey]['hist_shape'],
                                                                    dtype=object)
            for i in range(self.sysfcn[particle][systKey].shape[0]):
                for j in range(self.sysfcn[particle][systKey].shape[1]):
                    for k in range(self.sysfcn[particle][systKey].shape[2]):    
                        self.sysfcn[particle]['baseline_correction'][i,j,k] = \
                                                sf.dummy_func(align_sys_wsum[i,j,k]/align_sys_wcorrsum[i,j,k])


    def __init__(self, 
                 observables  = ['reco_energy', 'reco_zenith', 'pid'],
                 bin_edges   = [10**np.linspace(0.8,2.,11),
                                np.arccos(np.linspace(-1,0.,7))[::-1],
                                np.array([ -np.inf , 0.7, np.inf])],
                 user         = 'quickstart', 
                 LEaxis       = [],
                 weight_keys  = ['weight_e', 'weight_mu', 'weight_mu_pi'], 
                 break_energy = -1,
                 extra_cuts   = {},
                 table_nbins  = -1,
                 detsys_nuspecs = {'domeff': [1., 1],
                                   'hole_ice': [0.02, 2]},
                 detsys_muspecs = {'domeff': [1.0, 1]},
                 detsys_redo    = False,
                 detailed_detsys  = False,
                 verbose          = False,
                 use_kde_bg       = False,
                 use_kde_sys      = False,
                 quantile_mc      = (0., 1.),
                 default_settings = {},
                 sysfile          = 'DRAGON_detector_systematics.pckl'):


        # Store the arguments used
        self.iniDict = deepcopy(locals())
        self.iniDict['self']=None # module cannot be pickled
        del self.iniDict['self']
        print '\n ************************************ '
        print ' **** oscFit3D v1.1 - dataLoader ****' 
        print ' ************************************\n'
        print 'dataLoader: Initializing user', user
        
        # Importing user and settings
        self.user = __import__(user)
        reload(self.user)
        self.local_settings  = __import__('oscFit_settings')


        # Verify that the user has not repeated an observable
        if len(observables) > len(np.unique(observables)):
            print 'dataLoader: You are repeating an observable. You cannot do that. Exiting ...'
            sys.exit()
        # Check that the number of observables and edges match
        if len(observables) != len(bin_edges):
            print 'dataLoader: The number of observables and bin edges must match. Exiting ...'
            sys.exit()

        # Everything is done in 3D - If less observables are given, the code generates dummy dimensions
        # The hist_obs include these dummy dimensions
        num_observables = len(observables)
        if num_observables > 3:
            print 'dataLoader: For more than 3 observables use another sample (dataLoader instance)'
            sys.exit()
        elif num_observables < 3:
            print 'dataLoader: Dummy observable(s) created to use 3D histograms'
            hist_obs    = observables + [observables[0]]*(3-num_observables)
            bin_edges   += [np.array([-np.inf, np.inf])]*(3-num_observables)
        else:
            hist_obs    = observables

        # Declaration of variables used in dataLoader class
        self.sec2years   = 3600.*24*365
        self.bin_edges   = bin_edges
        self.observables = observables
        self.hist_obs    = hist_obs            # Includes dummy dimensions to keep everything 3D
        self.LEaxis          = LEaxis
        self.oscMC           = None
        self.we              = weight_keys[0]
        self.wmu             = weight_keys[1]
        self.wmu_pi          = weight_keys[2]
        self.verbose         = verbose
        self.ma_weights      = {}
        self.extra_cuts      = deepcopy(extra_cuts)
        self.break_energy    = break_energy
        self.detailed_detsys = detailed_detsys
        self.sysvar_params   = []
        self.detsys_nuspecs  = detsys_nuspecs
        self.detsys_muspecs  = detsys_muspecs
        self.unweighted_histograms = {}

        # Import KDE modules if they will be used at all
        if use_kde_bg or use_kde_sys:
            print 'dataLoader: KDE is set ON, importing KDE tools'
            from KDE_tools import smoothByKDE
            self.smoothByKDE = smoothByKDE

        # The limits on the observables (reco_energy, reco_zenith) are used as extra cuts - less computations
        # PID is not cut on (pid bias paramter). Using KDEs will override cuts on observables
        skip_observable_cut = ['new_pid', 'pid', 'pid2', 'LE']
        for i, one_observable in enumerate(self.observables):
            if one_observable in skip_observable_cut:
                continue
            self.extra_cuts[one_observable] = [np.min(self.bin_edges[i]), np.max(self.bin_edges[i])]

        # The L/E way needs this additional caveat - everything else should work as is. Case with/without PID.
        if len(self.LEaxis) != 0:
            if len(self.observables) == 3:
                self.observables = ['LE', self.observables[2]]
                self.hist_obs    =  self.observables + ['LE']
                self.bin_edges   = [ LEaxis, bin_edges[2],
                                    np.array([-np.inf, np.inf])]
            else:
                self.observables = ['LE']
                self.hist_obs    = ['LE']*3
                self.bin_edges = [ self.LEaxis] + [np.array([-np.inf, np.inf])]*2


        if table_nbins != None and table_nbins > 0:
            self.prepare_tables = True
        else:
            self.prepare_tables = False

        self.default_data_settings = oscfit_default_values.default_data_settings
        self.default_data_settings.update(  default_settings  ) 
        
        sysvar_skipkeys = ['norm_nu', 
                           'norm_noise', 
                           'atmmu_f',
                           'noise_f',
                           'domeff',
                           'hole_ice', 
                           'had_escale',
                           'atmmu_template', 
                           'norm_atmmu', 
                           'add_detector_systematics',
                           'correct_honda_flux']
        self.sysvar_keys = set(self.default_data_settings.keys()) - set(sysvar_skipkeys)

        self.oscCalc         = oscCalc.OscCalc(oscMode = self.default_data_settings['oscMode'],
                                               doTables = self.prepare_tables,
                                               nbins=table_nbins)
        self.nue   = {}
        self.numu  = {}
        self.nutau = {}
        self.atmmu = {}
        
        
        # Before loading all MC, check if there are systematic histograms for the settings requested
        baseline_only = False
        if os.path.exists(sysfile) and not detsys_redo:
            print 'dataLoader: Trying to load systematic histograms from file', sysfile
            infile = open(sysfile, 'r')
            self.detsys = pickle.load(infile)
            infile.close()
            baseline_only = True
        
        self.loadNeutrinoMC(baseline_only, use_kde_sys, quantile_mc)
        self.loadAtmMuMC(baseline_only, use_kde_bg, quantile_mc)

        if observables[0] != 'energy' and observables[1] != 'zenith':
            self.loadAtmMuData(use_kde=use_kde_bg)
            self.loadNoiseBackground(use_kde=use_kde_bg, 
                                     truth=self.default_data_settings,
                                     quantile_mc = quantile_mc)
            # Add a check of whether user *has* this value defined
            self.loadAtmMuDataAux(use_kde=use_kde_bg)
        
        # The Barlow likelihood needs unweighted histograms - noise calculated before. Get the rest.
        nuh = self.getNuHistogram(weight_power = 0.)
        self.unweighted_histograms.update(nuh)
        if bool(self.atmmu):  # Empty dictionaries evaluate to False
            self.unweighted_histograms['atmmu_histo'] = self.getAtmMuHistogram(weight_power=0.)

        # Get the histograms for different detector simulations
        self.nu_keys = ['nue','numu','nutau','nc']
        if not baseline_only:
            print 'dataLoader: Generating new systematic histograms'
            if use_kde_sys:
                print 'dataLoader: Using KDE tool for detector systematics. This will take a while ...'
            self.stageSystematicVariations()
            self.loadSystematicVariations(params = self.default_data_settings, 
                                          use_kde=use_kde_sys)

            print 'dataLoader: Done loading systematic variations - will dump to file now!'
            print 'File: ', sysfile
            pickle.dump(self.detsys, open(sysfile, 'w'))
            self.oscCalc = None
            reload(oscfit_default_values)

        # Build and combine the detector systematic functions
        self.buildSystematicFunctions()
        self.combineSystematics()
        self.default_data_settings = oscfit_default_values.default_data_settings
        self.default_data_settings.update(  default_settings  ) 
        #self.oscCalc = None
