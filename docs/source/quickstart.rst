.. _quickstart:

======================
 Quick start
======================


Anybody using oscFit should take some time to read at least the essential sections of this documentation. However, if you are eager to see it at work or want to test your installation, the following set of instructions will get you started:

1. **Obtain the latest stable version of oscFit**

   Go to :ref:`installation` and follow the instructions to get and configure the code.

2. **Get simulation files in oscFit-readable format**

   oscFit expects the user to provide simulated events in a very specific format. The instructions to do so from scratch can be found in :ref:`pickling`. For this example we will skip producing the files, and instead will use a previously converted simulation set. The set is stored in Madison::

     /data/user/jpyanez/oscFit/oscFit_quickstart_data.tar.gz

  Copy the file, place it in a suitable location and decompress it::

    tar -zxvf oscFit_quickstart_data.tar.gz

3. **Define a user**

   * Go to the */user_settings/* directory and open the *quickstart.py* file.

   * Find the *data_dir* variable, which should be the first one in the file, and write down the directory where you stored the simulation you just downloaded. Leave everything else as it is.

   * Save the changes
     

4. **Run the fit**

   * Go to */scripts/single_fits/*

   * If you *did not* install iminuit, uncomment the line that reads::
       
       #'minimizer':      'SLSQP',

     This will let the code use a scipy minimizer, which most users have installed by default.

   * Run::

       $ python quickstart_fit.py


The script will load the simulation, produce a set of pseudo-data and proceed to fit these data. Every minimization step will be printed out in your screen. Once the code finishes, the best fit values will be displayed. The whole process should take less than 5 minutes in a reasonable machine.

*Once you see a list of parameters and best fit values, you are done!*
