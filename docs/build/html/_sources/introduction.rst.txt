Introduction
============

What is pySiPM?
***************
pySiPM is a toolkit developed with the aim of simulating SiPM signals.
Even thou it has been initially developed to accomplish a specific task it has evolved in a way such that it can simulate generic SiPM signals.

Why it was created?
*******************
Originally it has been developed to simulate the SiPM digitized signals using the information given by a Geant4 simulation of the IDEA Dual Readout Calorimeter in which about 130 millions of SiPMs are coupled with scintillating and cherenkov fibers.
Its aim is to study the effects induced by the SiPM readout system in the detector response and eventually to help in the right choice of sensors and electronic setup.
It is also used to simulate the full geometry of IDEA Calorimeter and reconstruct high energy physics events with the granularity offered by the SiPM readout system.

Its aim is to help in the choice of the right SiPM sensor for a specific task and show how the stochastic noise effects may affect the values measured. Its aim is not to describe a SiPM in the most detailed way but to show how different parameters of the SiPM, signal shape and electric readout may affect the value measured by the user.


Who has created it?
*******************
| pySiPM has been developed for the IDEA collaboration at University of Insubria, department of Science and HighTech.
| The group that has developed pySiPM has an expertise in Silicon PhotoMultipliers devices and their applications.

* Edoardo Proserpio, physics student at Insubria University `eproserpio@studenti.uninsubria.it <eproserpio@studenti.uninsubria.it>`_
* PhD Massimiliano Antonello, research fellow at INFN Milano `massimiliano.antonello@mi.infn.it <massimiliano.antonello@mi.infn.it>`_
* PhD Romualdo Santoro, researcher and professor at Insubria University `romualdo.santoro@uninsubria.it <romualdo.santoro@uninsubria.it>`_
