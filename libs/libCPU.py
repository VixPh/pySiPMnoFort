# In this file I define all the functions I will use in the main file of simulation
from variables import *
from libs.FortranFunctions import signalgenfortran
###############################################################################################################
################>>>   EDITING THIS FILE MAY SERIOUSLY COMPROMISE SIMULATION BEHAVIOUR   <<<######################
#################################################################################################################

def PulseCPU(t,h):
	"""
	PulseCPU(t,h)

	Function that generates the signal from a single SiPM cell. This is the "full" version that computes the signal shape on CPU by evaluating the signal shape function.

	Parameters
	----------
	t : int32
		Time at which the cell is triggered
	h : float32
		The relative pulse height of the cell signal

	Returns
	-------
	s : np.ndarray
		Array containing the generated cell signal
	"""
	gainvar = np.float32(random.gauss(1,ccgv))											# Generate random ccgv for each fired cell
	s = signalgenfortran(t,h,tfall,trise,sigpts,gainvar)								# Calculate signal
	nap = poisson(ap)																	# Generate number of afterpulses
	if nap>0:																			# If there are afterpulses generate theyr signals like the pe ones
		for i in range(nap):
			apdel = random.expovariate(1/tauapfast)+random.expovariate(1/tauapslow)		# APs have a time delay exponential distribution
			tap = np.int32(apdel/sampling+t)											# Each afterpulse has a delay from its "main" signal that follows a exp distribution
			hap = 1-exp(-apdel/tfall)													# AP signal height depends exponentially by the delay
			s += signalgenfortran(tap,hap,tfall,trise,sigpts,gainvar)
	return(s)

def SiPMSignalAction(times,sigH,SNR,basespread):									# Function that passes signals times and height to main function for generating signals
	"""
	SiPMSignalAction(times,sigH,SNR,basespread)

	Function that passes signal height and times to the main function that generates single signals. Also adds noise.

	Parameters
	----------
	times : np.ndarray(int32)
		Array containing the time at wich SiPM cells are fired, including xt events (sorted)
	sigH : np.ndarray(float32)
		Array containing the pulse height of each fired SiPM cell
	SNR : double
		The signal to noise ratio of the noise to add
	basespread : double
		Sigma of the value to add as baseline

	Returns
	--------
	signal : np.ndarray
		Array containing the generated SiPM signal
	"""
	baseline = random.gauss(0,basespread)					# Add a baseline to the signal (gaussian distribution)
	signal = np.random.normal(baseline,SNR,sigpts)			# Start with gaussian noise
	for i in range(times.size):
		signal += PulseCPU(times[i],sigH[i])				# Generate signals
	return(signal)
