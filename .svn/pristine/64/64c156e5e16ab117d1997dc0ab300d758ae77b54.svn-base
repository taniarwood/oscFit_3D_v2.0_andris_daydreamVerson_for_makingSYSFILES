
==========================
 Loading data (dataLoader)
==========================

The dataLoader module in oscFit does all the heavy lifting. It will follow the instructions given by the user to load all the necessary information (such as simulation files) in order to produce histograms of how the data should look like given a set of parameters. Most of the code in this module is dedicated to re-weighting simulation and interpreting different sets to predict the detector behavior for values that have not been simulated, like arbitrary DOM efficiencies. Users can take advantage of this module to understand a sample, even if a fit is not being performed.



Essentials: Basic usage
=======================

To import the module, do::

  import dataLoader

To use dataLoader, you will need to create an instance of it. An example of how you define an instance follows::

    loader =  dataLoader.dataLoader(
        user = 'newbie',
        observables =  ['reco_energy', 'reco_zenith', 'pid'],
        bin_edges   =  [10**np.linspace(0.8,1.75,9),
                        np.arccos(np.linspace(-1,0.,9))[::-1],
                        np.array([-np.inf, 0.7, np.inf])],
        weight_keys = ['weight_e', 'weight_mu']
	)

A list containing the meaning of all the arguments taken can be found in the advanced section.

If this is the first time you are declaring an instance with the arguments you have passed, dataLoader will proceed to load all the files defined in the user files. Once the simulation is loaded, the module will build the detection systematic functions. This is done by creating histograms of all the sets declared, and fitting a function that describes the impact of a given systematic uncertainty in each bin. 

Loading all the files and producing the histograms (reweighting each event in each set) can be pretty slow, taking up to a few minutes for large samples. To avoid this delay, the output histograms are stored so that they can be used next time you start a dataLoader instance. This means that, as long as you use the same arguments in dataLoader, the code will read the previously generated histograms instead of all the simulation samples. In this case, only the baseline simulation is loaded.


Producing pseudo-data
.....................

One of the most useful methods in dataLoader, from the user point of view, is 

.. currentmodule:: dataLoader
.. autoclass:: dataLoader()
   :members: loadMCasData

With this method, it is possible to produce pseudo-data with different parameters and understand their impact on the histogram that is analyzed. The output is a histogram that follows the binning giving when the instance of dataLoader was initialized. To get the pseudo-data::

  pseudo_data = loader.loadMCasData(my_settings)


The data settings
.................

The module *oscfit_default_values* contains a dictionary with all the values required by dataLoader to produce a histogram. A user can import the module and use the dictionary to create a prediction of how the data would look like under certain circumstances. For example, to scale the neutrino expectation to 10 years of data with 5% muon contamination::

  import oscfit_default_values

  my_settings = oscfit_default_values.default_data_settings
  my_settings['norm_nu'] = 10.
  my_settings['atmmu_f'] = 0.05

  pseudo_data = loader.loadMCasData(my_settings)

There are many parameters which affect the simulation prediction. The following table lists all of them, together with their default values and a brief explanation of what they do:

========================= ============== ================================================
Parameter                 Default value  Description
========================= ============== ================================================
dm31                      0.0025         Std. oscillations - 3-1 mass splitting (in eV^2) 
theta23                   0.74           Std. oscillations - 2-3 mixing angle (in radians)
theta13                   0.155          Std. oscillations - 1-3 mixing angle (in radians)
mix_angle                 1.0            sin(2theta)^2 in 2-flavor approximation 
dm41                      0.0            Sterile neutrino oscillation parameters - 4-1 mass splitting (in eV^2)
theta24                   0.0            Sterile neutrino oscillation parameters - 2-4 mixing angle (in radians)
theta34                   0.0            Sterile neutrino oscillation parameters - 3-4 mixing angle (in radians)
norm_nu                   1.0            Normalization of neutrion sample (in years)
norm_e                    1.0            Scaling factor of the nue/numu ratio in flux prediction
norm_tau                  1.0            Normalization of the tau apperance due to oscillations - should be 1 in std. oscillations
norm_noise                1.0            Normalization of the pure noise simulation - use either norm or fraction, not both
noise_f                   None           Fraction of the total number of events that are pure noise - see comment above
norm_nc                   1.0            Scaling factor of neutral current cross section in GENIE
norm_atmmu                0.0            Normalization of the atmospheric muon template given - use either norm or fraction, not both
atmmu_f                   None           Fraction of the total number of events that are atmospheric mouns - see comment above
uphor_ratio               0.0            Change in the up-to-horizontal ratio of neutrinos predicted in the flux
nubar_ratio               0.0            Change in the nu/nubar ratio from the flux prediction
nu_nubar                  1.0            Helicity to blame for the change in nu/nubar ratio - continuous value from -1 (nubar) to 1 (nu)
gamma                     0.0            Change in the spectral index of the neutrino flux prediction, E^(gamma)
axm_qe                    0.0            Modify the cross section of quasi-elastic interactions by changing the axial mass (in units of sigma)
axm_res                   0.0            Modify the cross section of resonant interactions by changing the axial mass (in units of sigma)
pid_bias                  0.0            Introduce a linear bias on the particle identification parameter - zero is no bias
domeff                    1.0            Optical efficiency of the DOMs with respect to simulation - can be larger than 1.0
hole_ice                  0.02           Scattering coefficient in the refrozen column of ice around the DOM (in 1/cm)
hi_fwd                    0.0            Modification of the angular acceptance of the DOM in the very forward region (head-on collisions)
had_escale                1.0            Scale the light output of hadronic cascades
atmmu_template            'mc'           Identifier of the template being used for atmospheric muons - options are 'data','mc','corsika'
simulation                'baseline'     Identifier of the neutrino simulation used to produce the expectation
oscMode                   'Vacuum'       Method for calculating oscillation probabilities: TwoNeutrino, Vacuum, Prob3, NuCraft, GLOBeS
oscTables                 False          Set to True to calculate probabilities maps instead of event-wise
ma_variations             True           Use the axial mass variations - only set to False if your simulation doesn't have them
add_detector_systematics  True           Modify the expectation using the polynomial functions fit from multiple simulation sets
correct_honda_flux        False          Apply Honda's correction to the neutrino flux from using AMS-II data in their cosmic ray spectrum
========================= ============== ================================================

Take special care when using normalizations and fractions for the atmospheric muons and the noise. Define either a fraction or a normalization, but not both.


Everything is 3D
.................

The dataLoader module internally works in three dimensions. If you have given one or two, the code will create the necessary number of dummy dimensions to reach three (details in advanced section). Users will notice this when exploring the histograms produced by the module. For example, if you use *reco_energy* as your only observable, and you do::
  
  print pseudo_data.shape
  '(8, 1, 1)'

where the 8 appears because that was the number of bins declared in the example given above. The extra dimensions have a single bin.

Visualizing histograms
......................

Since the data is simply a python array, you can use any python tools to visualize it. Histograms, however, do not look very nice in python. The oscFit project comes with a few custom functions to improve the displaying of histograms, stored in the module **jp_mpl**. Perhaps the most useful one is 

.. currentmodule:: jp_mpl
.. autofunction:: unfilledBar

To use it::

  import jp_mpl as jplot

  jplot.unfilledBar( loader.bin_edges[0],  # Edges of the histogram
                     pseudo_data[:,0,0],   # The bin content
		     color = 'blue')









Advanced: Understanding the module
==================================

Add the full dataLoader init here, with all the freaking options

.. currentmodule:: dataLoader()
.. autoclass:: dataLoader
   :members:

...
