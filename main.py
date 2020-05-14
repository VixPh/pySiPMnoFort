#Function of simulation
from libs.lib import *

def SiPM(times,other):
	"""
	SiPM(times,other)

	Function that calls all the procedures defined in libs to generate a complete SiPM event.

	Parameters
	-----------
	times : np.ndarray
		This array contains the arriving time of each photon on the SiPM sensor surface. This array is the input of the simulation.

	other : tuple
		This variable may contain other informations about the event generated. It can be the event id, the arriving time inserted in the simulation or the real number of photons inserted in the simulation. This tuple will be copied as it is in the output.

	Returns
	----------
	integral : double
		The integral of the signal calculated in the integration gate

	peak : double
		The height of the signal in the integration gate

	tstart : double
		The time of arrival of the signal in ns defined as the first sample over the threshld of 1.5

	other : tuple
		The copy of the `other` variable given in the input

	signal : np.ndarray(float32)
		If the options -W is enabled the complete SiPM signal will be passed in the output. Otherwise this output is None
	"""
	###SIMULATION CORE###
	dcrTime = addDCR(dcr)									#Generate DCR events (times)
	if dcrTime.size:
		times = hstack((times,dcrTime))						# Concatenate DCR events with photoelectrons events
		sigTimes,sigH = SiPMEventAction(times,xt)			# Update list of times considering SiPM matrix occupancy and recovery times
	else:
		sigTimes,sigH = SiPMEventAction(times,xt)			# Update list of times considering SiPM matrix occupancy and recovery times

	signal = SiPMSignalAction(sigTimes,sigH,SNR,basespread)	# Generate digital signals

	signalInGate = signal[intstart:intstart+intgate]		# Select signal in the integration gate
	integral = signalInGate.sum()*sampling					# Calculating integrals
	peak = signalInGate.max()								# Calculating peaks
	tstart = np.argmax(signalInGate>1.5)*sampling			# Calculating starting times of signals
	tovert = np.count_nonzero(signalInGate>1.5)*sampling	# Calculating time over threshld
	tpeak = np.argmax(signalInGate)*sampling				# Calculating peaking time
	if args.Graphics:
		if not args.signal:
			dev = 'cpu-fast'
		elif args.device =='cpu':
			dev = 'cpu'
		elif args.device == 'gpu':
			if (times.size < cpuThreshold)|(times.size > gpuMax):
				dev = 'gpu(cpu)'
			else:
				dev = 'gpu'
		sigPlot(signal,sigTimes,dcrTime,dev)
	if args.wavedump is None:
		signal = None
	return(integral,peak,tstart,tovert,tpeak,other,np.float32(signal))
