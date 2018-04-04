.. _userdef:

==================
Declaring a *user*
==================

The *user* in oscFit is a bit of a misnomer. A better name to this would be a **sample**, but we will stick to user for now.

A *user* is a declaration of all the files (simulation and data) which correspond to a given event selection. This includes neutrino simulation, atmospheric muon background simulation and/or data, pure noise simulation, and the experimental data that you want to fit to. The association between file name and the way these events were simulated is done in here. You can declare multiple users, for example, to test different event selections.


Essentials: create your user
============================

User files are stored in the */user_settigs* directory. Users are identified by their file names. To create your own user make a copy of the *user_example.py* file. Let's assume you want this user to be called *newbie*::

  cp user_example.py newbie.py

Now open the *newbie.py* file. You will need to modify the following fields:

* **Directory where all pickle files are stored**::

     data_dir = '/path/to/your/data/dir'

* **Neutrino simulation sets**

  The user definition assumes that neutrino simulation files have the name convention::

    <PDG code><RUN ID>.pckl (e.g. 14550.pckl)

  The run id's are used to populate the *mc_sets* dictionary. The first level of this dictionary defines the properties of the sets stored in it. At minimum, it should contain the *baseline* key (first line). Sets of systematic variations in the simulation, like DOM efficiency, should be stored under their corresponding identifying key. The content of that key is a dictionary that assigns simulation values to run id's. The snipet below exemplifies how this is done::

     mc_sets = { 'baseline' : {'baseline': '550'},
                 'domeff':    {<DOM eff value> : 'RUN_ID',
	                       <DOM eff value> : 'RUN ID', 
		    	       ...
			       },
                  'hole_ice':  {<Hole ice value>: 'RUN_ID',
	                       ...
			       },
	         ...
	         }

* **Atmospheric muon simulation sets**

  If you have atmospheric muon simulation sets, they should be declared following the same convention as for the neutrinos. A typical declaration looks like::

   atmmu_sets = {'baseline': {'baseline':'mymuons'},
                 # 'domeff': {1.04:'mymuons_104',
                 #            1.2:'mymuons_120',
                 #            1.3:'mymuons_130',
                 #            0.85:'mymuons_85',
                 #            1.:'mymuons_100'}
                 }

   Here the user has commented out the 'domeff' sets, so they will not be used.

  If you do not have any atmospheric muon simulation, declare an empty dictionary::

   atmmu_sets = {}


* **Atmospheric muon sets from data**

  In some cases your muon expectation is derived using data. The declaration of these files is much simpler::

     atmmu_data_files = [ '/path/to/file1.pckl', ...]

  If you do not have a muon expectation from data, declare an emtpy list::

     atmmu_data_files = []

* **Pure noise simulation**

  Pure noise simulation is declared as a list of files, in the same way as the muons from data above. Following the same convention, declare an empty list if you do not have pure noise simulation::
  
    pure_noise_files = ['/path/to/MyFileOfPureNoise.pckl']
    pure_noise_files = [] # No pure noise files


* **Experimental data** 

  Experimental data files are stored in a dictionary to split them by detector year::

    data = { '2011': '/path/to/2011_data.pckl',
             '2012': '/path/to/2012_data.pckl',
	      ...
	   }

  Leave the dictionary empty (data={}) in case you do not have any data.



Advanced: Mixing GENIE and NuGen
================================


.. todo::

   Write this section. Nobody in the low-energy group needs this at the moment, but it is still part of the code, so it's worth writing a paragraph or two.
