Setup
=====

External file setup
*******************
The recommended way to change the parameters of the simulation or the SiPM is to launch the simulation with the option :option:`-f filename.txt`. The file loaded will overwrite the default options.
You can decide to overwrite all the default options or just some of them.

.. note::

	The file loaded with :option:`-f filename.txt` is read as a .py file, so comments must be preceded by `#` and every line will be executed as Python code.

Here is an example of file for a Hamamatzu SiPM S13615-1025

.. literalinclude:: ../../files/HAM-S13615-1025.txt

| The list of editable variables can be found at ... with a detailed description at ...
| Configuration files for some SiPM models are present in the `files` folder and it is a good practice to put your configuration files here.

Manual setup
************
It is also possible to manually change the values in the `variables.py` file. Keep in mind that changes in this files are permanent.

.. danger::

	Editing the `variables.py` file may change the behavior of the simulation in an unexpected way! Edit this file only if you know what you are doing.

.. _clineopts:

Command line options
********************
Some command line options can be parsed to de simulation and can be used to quickly change some parameters of the simulation.

    -H  Shows an help message and exits. Using :option:`-h` will display the ROOT help message.
    -V  Shows the software version and exits.
    -q  Switch to silent operation of the simulation.
    -G  Displays each generated signal (only for debugging purposes).
    -g  Shows histograms of integral, peak and starting time at the end of the simulation.
    -NDCR  Turns off dark count events generation.
    -NXT  Turns off the generation of optical crosstalk events.
    -NAP  Turns off the generation of afterpulses events.
    -FDCR  Enables faster generation of dark counts :ref:`dcrtheory`
    -FXT  Enables faster generation of optical crosstalk :ref:`xttheory`
    -SIG  Enables computation of each single signal (slower)
    -j <N>  Set the number of parallel workers to N (default is the number of cores detected)
    -f <file.txt>  Reads the configuration file selected
    -w <file.root>  Writes the results on a rootfile.
    -W <file>  Saves each digitized wave in a group of a hdf5 file.
