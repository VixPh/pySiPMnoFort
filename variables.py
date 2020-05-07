# In this file I define all the global variables that I will use in other files.
# f = open('files/banner.txt')
# exec(f.read())

import importlib										# This module helps me to decide wich modules import and if they're installed
import numpy as np										# Numpy is used to handle arrays
if importlib.find_loader('cupy') is not None:
	import cupy as cp									# Module equivalent to numpy but on GPU (needs Nvidia CUDA to be installed https://developer.nvidia.com/cuda-downloads )
else:
	print('Cupy not found, unable to run on GPU, will use CPU')
import matplotlib.pyplot as plt							# Only for showing signals plots
import time, random, os, sys, argparse, struct, multiprocessing, matplotlib, uproot,ROOT,h5py	# Usefull modules for I/O and system configuration
from multiprocessing import Pool										# Module for multiprocessing routines
from numpy import exp, hstack, unique, sort, sum, cumsum				# Import Numpy functions to access them directly later (faster than calling numpy.submodule every time)
from numpy.random import randint, poisson, choice, exponential			# Random generation routines from Numpy
from matplotlib.ticker import MaxNLocator

plt.style.use('fast')													# Faster rendering of signals
plt.rc('lines',antialiased=False)
#################################################################################################################
##################################>>>  VARIABLES DEFINITIONS  <<<################################################
#################################################################################################################
global siglen		# Length of signal in ns
global sampling		# Samlping time in ns
global sigpts		# Total number of points in signal
global ncell		# Total number of cells in SiPM matrix
global cellside		# Number of cells in the side of SiPM matrix
global dcr			# Dark Count Rate in Hz
global xt			# Optical Cross Talk probability in %
global tfall		# Falling time of SiPM signal in ns
global trise		# Rising time of SiPM signal in ns
global cellrecovery	# Cell recovery time in ns
global intstart		# Start of integration gate in ns
global intgate		# Integration gate lenght in ns
global preg			# Lenght of pre-gate in ns
global SNR			# Signal to noise ratio in dB
global basespread	# Baseline spread (sigma)
global ccgv			# Cell to cell gain variation (sigma)
global normpe		# Normalization of peack height (1 pe  => peak  =  1)
global ap			# After pulsing probability in %
global tauapfast	# After pulses time distribution decay (fast) in ns
global tauapslow	# After pulses time distribution decay (slow) in ns
global fastDCR		# Enable fast generation of DCR (less precision)
global fastXT		# Enable fast generation of XT (less precision)
global cpuThreshold	# If there are more pe than this value swich to GPU
global gpuMax		# If there are more pe than this value swich back to CPU (else it exeeds memory of GPU)
#################################################################################################################
##################################>>>   EDITABLE VARIABLES   <<<#################################################
#################################################################################################################

###VARIABLES INITIALIZATION##

# Simulation parameters
siglen = 750			# in ns
sampling = 1			# in ns

# SiPM parameters
size = 1				# in mm
cellsize = 25			# in um
dcr = 200e3				# in Hz
xt = 0.03				# in %
ap = 0.02				# in %
tfall = 50				# in ns
trise = 1				# in ns
cellrecovery = 20		# in ns
tauapfast = 15			# in ns
tauapslow = 85			# in ns
SNR = 30				# in dB
basespread = 0.00		# relative to single pe peack height (sigma)	Off at the moment
ccgv = 0.05				# relative to single pe peack height (sigma)

# Trigger parameters
intstart = 20			# in ns
intgate = 300			# in ns
preg = 20			# in ns

# Simulation parameters
fastDCR = False
fastXT = False
fastSIG = True
cpuThreshold = 100
gpuMax = 2000

#################################################################################################################
#####################################>>> ARGUMENTSS PARSER  <<<##################################################
#################################################################################################################
description = '''Software for SiPM simulation'''
epilog = '''Try -g -G to get started :)\
\n\nDeveloped for IDEA Dual Readout Calorimeter collaboration\
\nEdoardo Proserpio:\teproserpio@studenti.uninsubria.it\
\nMassimiliano Antonello:\tmassimiliano.antonello@mi.infn.it\
\nRomualdo Santoro:\tromualdo.santoro@uninsubria.it'''
parser = argparse.ArgumentParser('pySiPM',add_help=False,description=description,
epilog=epilog,
formatter_class=argparse.RawDescriptionHelpFormatter)

# parser._optionals.title = 'Options for the simulation'
# parser.add_argument('-H','--help',action='help',default=argparse.SUPPRESS)
# parser.add_argument('-V', '--version', action='version',version='%(prog)s 0.1')
# parser.add_argument('-d','--device',nargs='?',type=str,help='Select device for signal generation',choices=['cpu','gpu'],default='cpu')
# parser.add_argument('-g','--graphics',action = 'count',help = 'Histograms of generated events')
# parser.add_argument('-q','--quiet',action = 'count',help = 'Quiet')
# parser.add_argument('-w','--write',nargs = '?',type = str,help = 'File to write as output',metavar='filename.root')
# parser.add_argument('-G','--Graphics',action = 'count',help = 'Plot each signal (For debug purposes only)')
# parser.add_argument('-j','--jobs',type = int,help = 'Number of jobs for multiprocessing',metavar='N')
# parser.add_argument('-NDCR','--nodcr',action = 'count',help = 'Set DCR rate to 0')
# parser.add_argument('-FDCR','--fastdcr',action = 'count',help = 'Faster generation of DCR')
# parser.add_argument('-NXT','--noxt',action = 'count',help = 'Set XT rate to 0')
# parser.add_argument('-FXT','--fastxt',action = 'count',help = 'Faster generation of XT')
# parser.add_argument('-NAP','--noap',action = 'count',help = 'Set AP rate to 0')
# parser.add_argument('-SIG','--signal',action='count',help = 'Generate each signal independently (slower)')
# parser.add_argument('-f','--fname',nargs = '?',type = str ,help = 'Configuration file',metavar='filename.txt')
# parser.add_argument('-W','--wavedump',nargs = '?',type = str ,help = 'Output Digitized Waveforms on hdf5 file',metavar='groupname')
# parser.add_argument('-D','--clean',action='count',help = 'Clear old output files')
# args = parser.parse_args()
# del epilog
# del description
#
# if importlib.find_loader('cupy') is None:
# 	args.device = 'cpu'						# If Cupy is not installed use CPU
#
# if args.jobs is not None:					# If a number of jobs is choosen then it is used
# 	nJobs = args.jobs
# else:
# 	nJobs = multiprocessing.cpu_count()		# If not specified all cores are used
#
# if args.quiet:
# 	sys.stdout  =  open('/dev/null','w')	# Move all printed output to inexisting file # TODO: Find a smarter way to do so
#
# if args.nodcr:								# Set dcr rate to 0
# 	dcr = 0
# 	fastDCR = True
#
# if args.noxt:								# Set xt rate to 0
# 	xt = 0
# 	fastXT = True
#
# if args.noap:								# Set ap rate to 0
# 	ap = 0
#
# if args.fastdcr:							# Sligtly faster generation of DCR
# 	fastDCR  =  True
#
# if args.fastxt:								# Sligtly faster generation of XT
# 	fastXT  =  True
#
# if args.signal:								# Refer to libs/lib file for a detailed description
# 	fastSIG = False
#
#
# if args.fname is not None:
# 	print(f'\nReding SiPM setting from: {args.fname:s}')
# 	print('___________________________________')
# 	f = open(args.fname)
# 	print(f.read()+'___________________________________\n')
# 	f.close()
# 	f = open(args.fname)
# 	exec(f.read())							# Eventually load setting from external file
# else:
# 	print('\nUsing default SiPM settings!\n')	# Use settings defined above
#
# print('Detected %d cores...\r'%(multiprocessing.cpu_count()))
# print('Initializing simulation on %d cores...\n'%(nJobs))
# if args.signal is None:
# 	print('Generating signals with the fast method on CPU (default)...')
#
# if args.signal is not None:
# 	if args.device is None:
# 		print('Generating signals on CPU')
# 	if args.device is not None:
# 		print('Generating signals on '+args.device.upper())
#
# if args.clean is not None:
# 	if os.path.exists('waveforms.hdf5'):
# 		os.remove('waveforms.hdf5')


#################################################################################################################
#################################>>> NOT EDITABLE VARIABLES   <<<################################################
#################################################################################################################
# Conversion of time units from ns to units of sampling times. All the simulation works using samples as time unit.
sampling *= 1.
tfall = np.float32(tfall/sampling)
trise = np.float32(trise/sampling)
cellrecovery = np.float32(cellrecovery/sampling)
tauapfast = np.float32(tauapfast/sampling)
tauapslow = np.float32(tauapslow/sampling)
sigpts = int(siglen/sampling)
cellside = int(size/(cellsize*1e-3))
ncell = int(cellside**2)-1
intstart = np.min(intstart,0)
if intgate+intstart>siglen:
	ingate = siglen-intstart
intstart = int(intstart/sampling)
intgate = int(intgate/sampling)
preg = preg/sampling
preg = int(intstart-preg)
b = tfall/trise
normpe = 1 #(b**(1/(1-b))-b**(1/((1/b)-1)))**-1				# Uncommenting this will force the signal height to be the same as the number of pe
SNR = np.sqrt(10**(-SNR/20))
