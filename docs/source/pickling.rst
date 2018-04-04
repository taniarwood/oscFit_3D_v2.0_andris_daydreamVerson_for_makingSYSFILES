.. _pickling:

===========================================
 Producing oscFit-readable simulation files
=========================================== 

OscFit expects the user to provide the simulation and data files in the form of python dictionaries which contain numpy arrays. The *pickle* module is used to store and read the files.  This decision means that users have to convert their I3 or ROOT files for them to be used. While this intermediate step can be annoying, the result is that only the necessary information is loaded onto memory, and that loading the files themselves can be done in a matter of seconds.


Essentials: Converting I3 files
===============================

.. todo:: Test that the file actually works

The */resources* directory provides a scripts that transforms I3 files into oscFit format *(pckl)*. The script needs to be run from within an icetray environment, and uses the following modules:

* genie_icetray

* NewNuFlux (or another atmospheric neutrino flux calculator)

The script assumes that:

* You have two observables: energy and zenith angle.

* The output of the reconstructions are stored in two I3Particles. One for the muon track, and one for the hadronic cascade at the interaction vertex.

* The direction (zenith angle) is taken from the muon track. The energy is the sum of track and cascade energy.

* The events inside a file are 70% neutrinos and 30% antineutrinos.

In order to run the script, follow these steps:

1. From an icetray environment try::

     from icecube import NewNuFlux

     from icecube import genie_icetray
     
   If the modules load, you can go forward.

2. Open the *pickle_from_I3.py* file in */resources* and follow steps 1-4:

   * Let the script know how your observables are called

   * Set the directory where your data and MC is located

   * Set the directory where your oscFit data and MC files will be stored

   * Estimate the maximum number of events per file

   * Rename the file as *my_pickle_from_I3.py* or something similar

   If you don't rename the file you might end up committing these changes.

3. Run the script::
     
     $ python my_pickle_from_I3.py <RUN_NR>

   where RUN_NR are the run ids of your MC and data (e.g. 14500, 12500, IC86_1). You have to run it once for each set that you want to convert.


Advanced: The oscFit format
===========================

For simulation, each file is expected to include *the entire simulated set of a given particle* at the final level of the event selection. This means, for example, that a single file should include all muon neutrinos simulated with the same set of parameters. The same goes for, say, electron neutrinos. Sets of the same particles but simulated with different detector/medium properties should be stored in independent files.

Each file contains a dictionary. The keys of the dictionary are the different properties available. Each key points to a numpy array where these properties are stored in order. This means that the *n*-th element of all arrays inside a file correspond to the same event.

The following keys should be available in simulation files:

========================= =============================================================================================================
Key                                                                                                                  Value
========================= =============================================================================================================	
*ptype*                   Neutrino PDG code (-12, 12, -14, 14, -16, 16)

*energy*                  True neutrino energy (GeV)

*zenith*                  True neutrino zenith angle (radians)

*interaction*             Type of interaction: CC = 1, NC = 2

*weight_e*                Contribution of the event with this energy and zenith angle to the electron neutrino flux

*weight_mu*               Contribution of the event with this energy and zenith angle to the muon neutrino flux

*scattering*              Scattering process: DIS=1, COH, RES, QEL

*ma_res*                  Coefficients of 2nd order fit to the event weight resulting from varying the axial mass (RES)

*ma_qe*                   Coefficients of 2nd order fit to the event weight resulting from varying the axial mass (RES)
========================= =============================================================================================================	



Both data and simulation files should have keys that refer to the observables that will be used in the analysis. The typical ones are:

* *reco_energy*

* *reco_zenith*

The user is free to use as many as wanted and to name them at will. A typical example is to include a particle identification parameter to discern between track-like and cascade-like events.
