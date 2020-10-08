# In this file I define all the global variables that I will use in other files.
import files.banner
import argparse
import importlib
import multiprocessing
import os
import random
import struct
import sys
import time
import warnings
from multiprocessing import Pool

import h5py
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.ticker import MaxNLocator
from numpy import cumsum, exp, hstack, sort, unique
from numpy.random import choice, exponential, poisson, randint

if importlib.util.find_spec('cupy'):
    import cupy as cp
else:
    print('Cupy not found, unable to run on GPU, will use CPU')

matplotlib.use('Qt5Agg')
plt.style.use('fast')
plt.rc('lines', antialiased=False)

warnings.simplefilter('always', DeprecationWarning)
warnings.simplefilter('always', ImportWarning)


###############################################################################
##################################>>>  VARIABLES DEFINITIONS  <<<##############
###############################################################################
global SIGLEN	   	# Length of signal in ns
global SAMPLING		# Samlping time in ns
global SIGPTS		# Total number of points in signal
global NCELL		# Total number of cells in SiPM matrix
global CELLSIDE		# Number of cells in the side of SiPM matrix
global DCR			# Dark Count Rate in Hz
global XT			# Optical Cross Talk probability in %
global TFALL		# Falling time of SiPM signal in ns
global TRISE		# Rising time of SiPM signal in ns
global CELLRECOVERY # Cell recovery time in ns
global INTSTART		# Start of integration gate in ns
global INTGATE		# Integration gate lenght in ns
global PREG			# Lenght of pre-gate in ns
global SNR			# Signal to noise ratio in dB
global BASESPREAD   # Baseline spread (sigma)
global CCGV			# Cell to cell gain variation (sigma)
global NORMPE		# Normalization of peack height (1 pe  => peak  =  1)
global AP			# After pulsing probability in %
global TAUAPFAST    # After pulses time distribution decay (fast) in ns
global TAUAPSLOW    # After pulses time distribution decay (slow) in ns
global FASTDCR		# Enable fast generation of DCR (less precision)
global FASTXT		# Enable fast generation of XT (less precision)
global CPUTHRESHOLD # If there are more pe than this value swich to GPU
global GPUMAX		# If there are more pe than this value swich back to CPU


###############################################################################
##################################>>>   EDITABLE VARIABLES   <<<###############
###############################################################################
# Signal parameters
SIGLEN = 500        # in ns
SAMPLING = 0.1      # in ns

# SiPM parameters
SIZE = 1			# in mm
CELLSIZE = 25		# in um
DCR = 200e3			# in Hz
XT = 0.01			# in %
AP = 0.03			# in %
TFALL = 50			# in ns
TRISE = 1			# in ns
CELLRECOVERY = 20   # in ns
TAUAPFAST = 15		# in ns
TAUAPSLOW = 85		# in ns
SNR = 30			# in dB single pe peack height (sigma) Off at the moment
BASESPREAD = 0.00
CCGV = 0.05			# relative to single pe peack height (sigma)

# Trigger parameters
INTSTART = 5		# in ns
INTGATE = 300		# in ns
PREG = 0			# in ns

# Simulation parameters
FASTDCR = False
FASTXT = False
FASTSIG = True
CPUTHRESHOLD = 100
GPUMAX = 2000


###############################################################################
##########################>>> ARGUMENTSS PARSER  <<<###########################
###############################################################################
description = '''Software for SiPM simulation'''
epilog = '''Try -g -G to get started :)\
\n\nDeveloped for IDEA Dual Readout Calorimeter collaboration\
\nEdoardo Proserpio:\teproserpio@studenti.uninsubria.it\
\nMassimiliano Antonello:\tmassimiliano.antonello@mi.infn.it\
\nRomualdo Santoro:\tromualdo.santoro@uninsubria.it'''
parser = argparse.ArgumentParser('pySiPM',
                                 add_help=False,
                                 description=description,
                                 epilog=epilog,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser._optionals.title = 'Options for the simulation'
parser.add_argument('-H', '--help', action='help', default=argparse.SUPPRESS)
parser.add_argument('-V', '--version', action='version',
                    version='%(prog)s 0.2')
parser.add_argument('-d', '--device', nargs='?', type=str,
                    help='Select device for signal generation',
                    choices=['cpu', 'gpu'], default='cpu')
parser.add_argument('-g', '--graphics', action='count',
                    help='Histograms of generated events')
parser.add_argument('-q', '--quiet', action='count', help='Quiet')
parser.add_argument('-w', '--write', nargs='?', type=str,
                    help='File to write as output', metavar='filename.root')
parser.add_argument('-G', '--Graphics', action='store',
                    help='Plot each signal (For debug purposes only) specify interval in ms')
parser.add_argument('-j', '--jobs', type=int,
                    help='Number of jobs for multiprocessing', metavar='N')
parser.add_argument('-NDCR', '--nodcr', action='count',
                    help='Set DCR rate to 0')
parser.add_argument('-FDCR', '--FASTDCR', action='count',
                    help='Faster generation of DCR')
parser.add_argument('-NXT', '--noxt', action='count', help='Set XT rate to 0')
parser.add_argument('-FXT', '--FASTXT', action='count',
                    help='Faster generation of XT')
parser.add_argument('-NAP', '--noap', action='count', help='Set AP rate to 0')
parser.add_argument('-SIG', '--signal', action='count',
                    help='Generate each signal independently (slower)')
parser.add_argument('-f', '--fname', nargs='?', type=str,
                    help='Configuration file', metavar='filename.txt')
parser.add_argument('-W', '--wavedump', nargs='?', type=str,
                    help='Output Digitized Waveforms on hdf5 file', metavar='groupname')
parser.add_argument('-D', '--clean', action='count',
                    help='Clear old output files')

args = parser.parse_args()

if importlib.util.find_spec('cupy') is None:
    args.device = 'cpu'

if args.jobs:
    if isinstance(args.jobs, int):
        nJobs = args.jobs
    else:
        warnings.warn(f'Invalid entry for number of jobs: {args.jobs}. Using only one', category=UserWarning)
        nJobs = 1
else:
    nJobs = multiprocessing.cpu_count()  # If not specified all cores are used

if args.quiet:
    sys.stdout = open('/dev/null', 'w')

if args.nodcr:		# Set DCR rate to 0
    DCR = 0
    FASTDCR = True

if args.noxt:		# Set XT rate to 0
    XT = 0
    FASTXT = True

if args.noap:		# Set AP rate to 0
    AP = 0

if args.FASTDCR:  # Sligtly faster generation of DCR
    FASTDCR = True

if args.FASTXT:		# Sligtly faster generation of XT
    FASTXT = True

if args.signal:		# Refer to libs/lib file for a detailed description
    FASTSIG = False

if args.fname:
    print(f'\nReding SiPM setting from: {args.fname:s}')
    print('___________________________________')
    f = open(args.fname)
    [print(l.rstrip()) for l in f.readlines() if not l.startswith('#')]
    print('___________________________________\n')
    f.close()
    f = open(args.fname)
    exec(f.read())			# Eventually load setting from external file
else:
    print('\nUsing default SiPM settings!\n')  # Use settings defined above

print('Detected %d cores...\r' % (multiprocessing.cpu_count()))
print('Initializing simulation on %d cores...\n' % (nJobs))
if args.signal is None:
    print('Generating signals with the fast method on CPU (default)...')

if args.signal:
    if args.device is None:
        print('Generating signals on CPU')
    if args.device:
        print('Generating signals on ' + args.device.upper())
        if args.device == 'gpu':
            warnings.warn('Signal generation on GPU is deprecated... use CPU preferably',category=DeprecationWarning,stacklevel=3)

if args.clean is not None:
    if os.path.exists('waveforms.hdf5'):
        os.remove('waveforms.hdf5')

###############################################################################
#######################>>> NOT EDITABLE VARIABLES   <<<########################
###############################################################################
# Conversion of time units from ns to units of sampling times.
SAMPLING *= 1.
TFALL = np.float32(TFALL / SAMPLING)
TRISE = np.float32(TRISE / SAMPLING)
CELLRECOVERY = np.float32(CELLRECOVERY / SAMPLING)
TAUAPFAST = np.float32(TAUAPFAST / SAMPLING)
TAUAPSLOW = np.float32(TAUAPSLOW / SAMPLING)
SIGPTS = int(SIGLEN / SAMPLING)
CELLSIDE = int(SIZE / (CELLSIZE * 1e-3))
NCELL = int(CELLSIDE**2) - 1
if INTGATE + INTSTART > SIGLEN:
    warnings.warn(f'Integration gate of {INTGATE:.0f} ns exeeds signal length of {SIGLEN:.0f} ns',
                  category=UserWarning)
    INTGATE = SIGLEN - INTSTART
INTSTART = int(INTSTART / SAMPLING)
INTGATE = int(INTGATE / SAMPLING)
PREG = PREG / SAMPLING
PREG = int(INTSTART - PREG)
b = TFALL / TRISE
NORMPE = 1  # (b**(1/(1-b))-b**(1/((1/b)-1)))**-1
SNR = np.sqrt(10**(-SNR / 20))
