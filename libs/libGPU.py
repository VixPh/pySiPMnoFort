# In this file I define all the functions I will use in the main file of simulation
from variables import *
from libs.FortranFunctions import signalgenfortran

#################################################################################################################
################>>>   EDITING THIS FILE MAY SERIOUSLY COMPROMISE SIMULATION BEHAVIOUR   <<<######################
#################################################################################################################

###GENERATION OF SIGNALS###
signalShapeGPU = cp.ElementwiseKernel(											# CUDA kernel that generates the signal shape (exp(-t/tf)-exp(t/tr))
'int32 x, float32 tfall, float32 trise, float32 ccgv, float32 h',				# All signals are generated in a matrix,each row is a signal and then summed up column-wise
'float32 z',
'z = h*ccgv*(__expf(-x/tfall)-__expf(-x/trise));',
'signalShape')

def PulseGPU(t,h):
	"""
	PulseCPU(t,h)

	Function that generates the signal from all SiPM cells at once. This is the "full" version that computes the signal shape on GPU

	Parameters
	----------
	t : np.ndarray(int32)
		Array containing times at which each the cell is triggered
	h : np.ndarray(float32)
		Array containing the relative pulse height of each cell signal

	Returns
	-------
	s : np.ndarray
		Array containing the generated SiPM signal
	"""
	n = t.size																	# Number of signals to generate
	nap = poisson(ap*n)															# Generate number of afterpulses: ap per each signal so it ends up to be n*ap
	vect = (cp.arange(sigpts,dtype='int32')+cp.zeros((n,1),dtype='int32')-t[:,None])	# Generate matrix containing times of each fired cell (x variable of exponential function)
	vect[vect<0] = 0															# Zero before the signal
	gainvar = cp.random.normal(1,ccgv,(n,1),dtype='float32')					# Generate random ccgv for each fired cell
	h = h[:,None].astype('float32')												# Transpose array of height values (each generated signal has it's height)
	s = normpe*cp.sum(signalShapeGPU(vect,tfall,trise,gainvar,h),axis=0)		# Call kernel to generate singal
	if nap>0:																	# If there are afterpulses generate theyr signals like the pe ones
		apdel = cp.random.exponential(tauapfast,nap,dtype='float32')+cp.random.exponential(tauapslow,nap,dtype='float32')				# APs have a time delay exponential distribution
		apSig = randint(0,n,dtype='int32')										# Select wich signals will have ap
		tap = (apdel/sampling).astype('int32')+t[apSig]							# Each afterpulse has a delay from its "main" signal that follows a exp distribution
		hap = 1-cp.exp(-apdel/tfall)											# The pulse height is calculated considering RC circuit charge profile
		hap = hap[:,None].astype('float32')
		gainvar = cp.random.normal(1,ccgv,(nap,1),dtype='float32')				# Generate gain variation for ap
		vect = (cp.arange(sigpts,dtype='int32')+cp.zeros((nap,1),dtype='int32')-tap[:,None])
		vect[vect<0] = 0
		s += normpe*cp.sum(signalShapeGPU(vect,tfall,trise,gainvar,hap),axis=0)	# Add each afterpulse signal to main signal
	return(s)

# Signals are made with the two exponentials approximation and each peack height is smeared by a gaussian distribution (Cell to cell gain variation), each signals is also mutiplied by its height considering cell recovery time (h=1 if i's the firs event or less if the cell is fired multiple times)
def PulseCPU(t,h):
	"""
	PulseCPU(t,h)

	Function that generates the signal from a single SiPM cell. This is the "full" version that computes the signal shape on CPU

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
	s = signalgenfortran(t,h,tfall,trise,sigpts,gainvar)
	nap = poisson(ap)																	# Generate number of afterpulses
	if nap>0:																			#If there are afterpulses generate theyr signals like the pe ones
		for i in range(nap):
			apdel = random.expovariate(1/tauapfast)+random.expovariate(1/tauapslow)		# APs have a time delay exponential distribution
			tap = np.int32(apdel/sampling+t)											# Each afterpulse has a delay from its "main" signal that follows a exp distribution
			hap = 1-exp(-apdel/tfall)													# AP signal height depends exponentially by the delay
			s += signalgenfortran(tap,hap,tfall,trise,sigpts,gainvar)
	return(s)

def signalGen(times,sigH,SNR,basespread):												# Function that passes signals times and height to main function for generating signals
	"""
	signalGen(times,sigH,SNR,basespread)

	Function that passes signal height and times to the main function that generates signals. Also adds noise.

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
	baseline = random.gauss(0,basespread)											# Add a baseline to the signal (gaussian distribution)
	if (times.size < cpuThreshold) or (times.size > gpuMax):							# If there is low light generate on cpu (if there is too much light it wont fit on VRAM so fall back on cpu)
		signal = np.random.normal(baseline,SNR,sigpts)								# Add white gaussian noise to the signal (gaussian distribution)
		for i in range(times.size):
			signal += PulseCPU(times[i],sigH[i])									# Generate signals one by one and sum them
		return(signal)
	else:
		signal = cp.random.normal(baseline,SNR,sigpts,dtype='float32')				# Add white gaussian noise to the signal (gaussian distribution)
		signal += PulseGPU(cp.asarray(times),cp.asarray(sigH))						# Generate signals all at once
		return(cp.asnumpy(signal))
