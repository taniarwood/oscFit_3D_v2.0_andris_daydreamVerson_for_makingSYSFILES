================
 Introduction
================


What is oscFit?
===========================

The oscFit project is a collection of tools whose main goal is to compare atmospheric neutrino data to simulation, and determine the parameters in simulation that result in the best fit to the data. While the task is conceptually simple, there are mutliple steps in between which have to be solved for this comparison to be made, and oscFit offers a way to deal with all (or most) of them. 

..
   Amongst the trivial tasks are loading the data and simulation to memory, or producing and plotting histograms. As the reader will see, that part is implemented in a straight-forward manner in oscFit. The most complex task is dealing with the calculation of event weights. Some of this can be counter-intuitive, and so special attention is given to that topic throughout the documentation.


oscFit at work
===============================

In oscFit data and simulation are compared using histograms in three dimensions, which are populated using the simulated events and the detected ones. By default, histograms are compared using a *Likelihood*, although other modes are also present.

The code assumes that the simulation is *weighted*, i.e. events do not contribute equally in the resulting histogram of the simulation. These weights are calculated internally in order to produce an expectation which is based on values different from those used to produce the simulation. In oscFit you can easily answer questions such as: *how would the data look like if the axial mass used to model coherent scattering would be off by -1 standard deviations?*.

In order to find the set of simulation parameters which resembles the data the most, oscFit uses a *fitter*. For historical reasons, the default one is an implementation of MINUIT [#]_ for python, called iminuit [#]_ . Native python minmizers have also been implemented. Users should optimize the fitter settings to help it reach an answer in a reasonable number of comparisons.

The code is split into *modules*. The modules and their tasks are listed below:

* dataLoader - Load simulation and data to memory. Reweight simulation.
* oscCalc - Interface between calculators of oscillation probabilities and oscFit
* oscFit - Fit the simulation to data using a minimizer
* plotTools - Plot results and comparisons between data and simulation

Some of the modules, particularly *dataLoader*, have functions that can make the user's life easier. Explore them before you engage in writing complicated code to manipulate the simulation.





.. [#] https://en.wikipedia.org/wiki/MINUIT
.. [#] https://pypi.python.org/pypi/iminuit
