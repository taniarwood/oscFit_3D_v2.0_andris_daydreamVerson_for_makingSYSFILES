.. _fitting:

==================
Making a fit
==================

A fit is done in three steps:

1. Starting an instance of dataLoader

2. Obtaining the (pseudo-)data that the MC will try to reproduce

3. Setting the fit parameters and running oscFit

While it all can be done in an interactive console, users are encouraged to write this down in a script. 

Essentials: A simple fit
========================

In order to do a fit, you need to have *oscFit-readable* files and a *user* which knows where the files are and how they should be interpreted. To get those, visit :ref:`pickling` and :ref:`userdef`. If you have them, then proceed to the directory */scripts/single_fits*. Make a copy of the example file::

  cp oscFit_example.py oscFit_newbie.py

Open the oscFit_newbie.py. The file contains the minimum set of instructions to make a fit. This will give you a basic understanding of how things work. Additional options are described in the *Advanced* section below.

The first important piece of code in the file is the declaration of a *dataLoader* instance, which loads the files to memory::

    loader =  dataLoader.dataLoader(
        user = 'newbie',
        observables =  ['reco_energy', 'reco_zenith', 'pid'],
        bin_edges   =  [10**np.linspace(0.8,1.75,9),
                        np.arccos(np.linspace(-1,0.,9))[::-1],
                        np.array([-np.inf, 0.7, np.inf])],
        weight_keys = ['weight_e', 'weight_mu']
	)

The basic arguments that dataLoader takes are:

=================== ==============================================================
Argument            Allowed values
=================== ==============================================================
user                Name of the user file that will be used (for our example this is "newbie")
observables         List with the name of your observables in the pickle files. No more than three.
bin_edges           List of bin edges to be used in the histogram for each observable, in the same order as given in the *observables* argument (bin edges should be given as numpy arrays).
weight_keys         List of the keys that contain the atmospheric fluxe times neutrino simulation weight for electron and muon neutrino fluxes, respectively.
=================== ==============================================================

The next step is to create the (pseudo-)data that the fitter will try to match. If you are using **real detector data**, then the instruction is rather simple::

  data_histogram = loader.loadData(year = 'all')

where the argument year is a string of the year to be analyzed, or 'all' to use all years declared.

If you are using **pseudo-data** from the simulation itself, then you might want to first define the parameters that govern such data. There is a set of default values, defined in the *oscfit_default_values.py* file. If you define a value in your script, the default value will be substituted by the one you have chosen. You can also chose to use only default values. The most common ones that you might want to modify are::

  data_settings = {'dm31':      0.0025,      # In eV^2
                   'theta23':   0.7,         # In radians
		   'theta13':   0.155,       # In radians
		   'norm_nu':   3.,          # Normalization of the sample, in years
		   'gamma':     0.03,        # Modification of the spectral index of the flux
		   'atmmu_template': 'mc',   # Source of your muon template (mc, corsika, data) 
		   'norm_atmmu':0.1,         # Scaling factor of the muon template
		   'simulation':'baseline',  # Simulation set from which the data is taken
		   'oscMode':   'Vacuum',    # Mode to calculate oscillation probabilities
		  }

Once you have defined your settings, retrieve a pseudo-data histogram with::

  data_histogram = loader.loadMCasData(data_settings,
                                       statistical_fluctuations=False)

If you set the *statistical_fluctuations* argument to False the histogram you obtain will be the perfect expectation, also known as the Asimov ensemble. If you set it to True you will obtain a histogram with appropiate poisson fluctuations in each bin.

The next thing that should be defined are the **fit settings**. The entire list of possible settings is found in *oscfit_default_values.py*. Just as before, default values are used for settings that are not defined by the user. The fit settings is a dictionary with the structure::

  fit_settings = { 'simulation':   'baseline',           # MC set to be used for the fit
                   'dm31':         [0.002, False, 'NH'], # Starting value, True to fix, NH/IH (hierarchy)
		   'theta23':      [0.7, False],         # Starting value, True to fix
		   ...
		   'detector_syst': True,                # Set to False to avoid using detector systematics
		   'include_priors':True,                # False to avoid using priors
		   'printMode':     -1,                  # See the evolution of the fit (-1, 0, 1, 2)
		 }

The only thing that remains is to run the fit. Start an instance of oscFit.fitOscParams(), and pass the data you want to fit, as well as the loader and your settings::

  fitter = oscFit.fitOscParams()
  result = fitter(data_histograms = [data_histogram],
                  data_loaders    = [loader],
		  fit_settings    = fit_settings)

And now, you have to sit back and wait for the fitter to find you a minimum.


Advanced: More fit options
==========================

Defining your priors
--------------------

The likelihood calculation in oscFit applies penalty factors to the nuisance parameters in the fit when they acquire a value which differs from their best known one. The penalty terms assume the error on the nuisance parameter is Gaussian, and are of the form

:math:`\dfrac{(\eta - \eta_\mathrm{fit})^2}{\sigma_\eta^2}`.

A set of values for the best known value and the error on the nuisance parameters (also called *priors*) are given in the code. The user can override these values by passing a dictionary with the values desired to the fitter. For example::

  my_own_priors = {'norm_e': [ 1., 0.05 ],
                   'gamma':  [ 0.05, 0.1],
		   'domeff': [ 1.04, 0.06]
		   }

  fitter = oscFit.fitOscParams()

  result = fitter(data_histograms = data_histograms,
                  data_loaders    = data_loaders,
		  fit_priors      = my_own_priors}


Store fit results
-----------------

The following code snippet will store the result of your fit in the directory defined in the oscFit_settings.py file in the user_settings directory::

  result_dict = {'result': result, 
                 'fit_settings': fit_settings,
		 'loader_settings': [x.iniDict for x in data_loaders],
		 'data_histograms': data_histograms,
		 'data_settings':   data_settings


.. todo:: I need to standardize this so the user cannot fuck up. I should save it as a possibility inside the fitter function. No need to store stuff afterwards.


Running with HESSE errors (MINUIT)
----------------------------------


Store fit settings
------------------

Fit with different MC
---------------------

L/E fits
--------


Including more observables
--------------------------

Fitting multiple samples
------------------------


Defining your priors
--------------------


Starting values
---------------


Legacy mode
-----------

Detector systematics
--------------------
