Files Structure
===============
In this section the file structure used to build the simulation is explained. You can find a brief description of the files used and they content. For a detailed explanation see the individual sections.


Main
-----

The simulation is made of a series of function and variables that are used in the main file to generate a SiPM event. In the main file the number of cells triggered is calculated and the signal is generated accordingly. In this file the signal is also analysed.

The structure of the main file may be modified by the user to include more advanced analysis of the digitised.  signal. 

main.py
*******
In this file the functions present in the `libs` folder are combined to generate a complete SiPM event

.. currentmodule:: main
.. autofunction:: SiPM

Variables
----------

variables.py
*************
| The file `variables.py` contains all the variables needed to describe the SiPM sensor, the simulation behavior and the signal analysis.
| Here is a brief description of all the parameters with their default value.

* SiPM Settings:

  * :code:`size = 1` Size of the SiPM sensitive area in mm (side)
  * :code:`cellsize = 25` Pitch of a single sensitive cell in um
  * :code:`dcr = 200e3` Dark count rate in kHz
  * :code:`xt = 0.03` Optical crosstalk probability
  * :code:`ap = 0.02` Afterpulsing probability
  * :code:`trise = 1` Rising time of the signal in ns
  * :code:`tfall = 50` Falling time of the signal in ns
  * :code:`cellrecovery = 20` Recovery time of the single SiPM cell in ns
  * :code:`tauapfast = 15` Time constant of afterpulses delay in ns (fast component)
  * :code:`tauapslow = 85` Time constant of afterpulses delay in ns (slow component)
  * :code:`ccgv = 0.05` Cell to cell gain variation (sigma)
  * :code:`SNR = 30` Signal to noise ratio
  * :code:`basespread = 0` Spread of the baseline of the signal (sigma)

* Signal analysis settings:

  * :code:`siglen = 500` Total length of the signal in ns
  * :code:`sampling = 0.1` Samlping time in ns (1 = 1ns, 0.1 = 100ps)
  * :code:`intstart = 20` Start of the integration gate in ns
  * :code:`intgate = 300` Length of the integration gate in ns
  * :code:`pregate = 20` Length of the pregate

* Simulation settings:

  * :code:`fastDCR = False` Enables faster generation of dark count events
  * :code:`fastXT = False`  Enable faster generation of crosstalk events
  * :code:`fastSIG = True`  Enables faster computation of signals

For a more detailed description of each variable see ...

Functions
----------
libs/lib.py
***********

| This file contains all the functions that are used to generate the event and the SiPM signal.
| Here is a brief description of all the functions with their signature.

.. currentmodule:: libs.lib
.. autofunction:: addDCR
.. autofunction:: evtsGen
.. autofunction:: signalGen
.. autofunction:: PulseCPU
.. autofunction:: sigPlot
.. autofunction:: initializeRandomPool
.. autofunction:: somestats

libs/libCPU.py
**************

| This file contains the functions used to compute the SiPM signal using the CPU.
| Here is a brief description of all the functions with their signature.

.. currentmodule:: libs.libCPU
.. autofunction:: signalGen
.. autofunction:: PulseCPU

libs/libGPU.py
**************

| This file contains the functions used to compute the SiPM signal using the GPU.
| Here is a brief description of all the functions with their signature.

.. currentmodule:: libs.libGPU
.. autofunction:: signalGen
.. autofunction:: PulseCPU
.. autofunction:: PulseGPU

Other
------

scr
***

This folder contains the source file of the FORTRAN functions that are used to speed-up the simulation

files
*****

This folder contains some files used by the simulation and contains pre-configurations of some SiPM models that can be eventually loaded.
For an explanation on how to write and load an external configuration file see ...
