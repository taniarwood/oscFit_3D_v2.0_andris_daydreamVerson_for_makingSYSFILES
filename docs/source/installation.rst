.. _installation:

================
 Installation
================


Pre-requisites
==============

The code was developed using

* Python 2.7.2
* NumPy 1.9.0
* SciPy 0.11.0


Minimizers
..........

The minimizer guides the comparisons between data and simulation, and decides when the best possible fit has been reached. There are multiple strategies for doing this, and their effectiveness depends on the problem at hand. In its default mode, oscFit uses the MIGRAD algorithm from the `iminuit <https://pypi.python.org/pypi/iminuit>`_ package (follow the link for installation instructions). Some routines implemented in the scipy package are also supported. They have not been tested to the same extent, but so far they seem promising.

To install iminuit using pip install::

  pip install cython
  pip install iminuit


Oscillation probabilities
.........................

Atmospheric neutrinos oscillations are affected by the Earth's matter density profile. A short explanation of the observable effects above a few GeV, in the context of neutrino telescopes, can be found in `this paper <http://arxiv.org/abs/1509.08404>`_. The calculation of oscillations with matter effects is non-trivial, and the code to do this numerically has been developed multiple times. This is why the oscFit project uses third party software for calculating oscillation probabilities.

The project can currently interface with three calculators that perform neutrino oscillations calculations with matter effects. In order to use them, these calculators need to be installed by the user. They are:

* `Prob3++ <http://www.phy.duke.edu/~raw22/public/Prob3++/>`_ -- The code is written in C++. Requires a python wrapper. Currently using the one from PISA.
* `GLoBES <https://www.mpi-hd.mpg.de/personalhomes/globes/>`_ -- The code is written in C++. Using a python wrapper written by A. Terliuk. Can handle sterile neutrinos.
* `NuCraft <https://www.mpi-hd.mpg.de/personalhomes/globes/>`_ -- Python code. Handles sterile neutrinos.

The calculators are listed by how fast they perform. If you are interested in standard neutrino oscillations (no sterile states), Prob3++ is the only one you will need.

Neutrino oscillations in vacuum do not need to be solved numerically. This is why oscFit already contains a fast method to calculate them, either in the two or three flavor mode. Because of it, you can start using oscFit without installing any of the packages described above. However, assuming that neutrinos travel through vacuum as they cross the Earth is not a very good approximation, so this method should not be used to compare with experimental data.


Obtaining and configuring the code
==================================

The oscFit code lives in this svn repository:
  http://code.icecube.wisc.edu/svn/sandbox/yanez/oscFit/
Stable versions of the code can be found under */releases*. The most up-to-date version lives in */trunk*.

After getting a version of the code:

0. If you have a previous version, just do::

     ./setup.py

   and follow the instructions. You are done! Otherwise, continue with this list.

1. Go to the *user_settings* folder and execute::

     cp oscFit_settings_template.py oscFit_settings.py


2. Open the copied file **oscFit_settings.py**

3. Set the directory where your fits (fits_dir) and scans (scans_dir) will be stored.

4. Set the directory where your auxiliary files will be stored (kde_database, detsys_database).

5. If you have installed external oscillations calculators, store their location in the variables provided. If you did not install any of the external calculators, leave them as they are. Save the changes.

6. Whenever you want to use any of the modules of oscFit, you will need to tell python where the project lives. Example scripts do it with these three lines at the beginning of the file::

     import os, sys
     modules_dir = os.path.dirname(os.path.realpath(__file__)).split('/scripts')[0] + '/modules'
     sys.path.append(modules_dir)

 If you use the directory structure provided by oscFit you can copy/paste these lines into your new scripts.

