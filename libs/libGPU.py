# In this file I define all the functions I will use in the main file of simulation
from libs.FortranFunctions import signalgenfortran
from variables import *

###############################################################################
##>>>   EDITING THIS FILE MAY SERIOUSLY COMPROMISE SIMULATION BEHAVIOUR   <<<##
###############################################################################

###GENERATION OF SIGNALS###
signalShapeGPU = cp.ElementwiseKernel(
    # CUDA kernel that generates the signal
    # Signals are generated in a matrix,each row is a signal, summed up column-wise
    'int32 x, float32 TFALL, float32 TRISE, float32 CCGV, float32 h',
    'float32 z',
    'z = h*CCGV*(__expf(-x/TFALL)-__expf(-x/TRISE));',
    'signalShape')


def PulseGPU(t, h):
    """
    PulseCPU(t,h)

    Function that generates the signal from all SiPM cells at once.
    This is the "full" version that computes the signal shape on GPU.

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
    n = t.size  # Number of signals to generate
    nap = poisson(AP * n)    # Number of afterpulses
    # Generate matrix containing times of each fired cell
    vect = (cp.arange(SIGPTS, dtype='int32') + cp.zeros((n, 1), dtype='int32') - t[:, None])
    vect[vect < 0] = 0   # Zero before the signal
    # Generate random ccgv
    gainvar = cp.random.normal(1, CCGV, (n, 1), dtype='float32')
    h = h[:, None].astype('float32')    # Transpose array of height values
    # Call kernel to generate singal
    sig = normpe * cp.sum(signalShapeGPU(vect, TFALL, TRISE, gainvar, h), axis=0)
    # If there are afterpulses generate theyr signals
    if nap:
        apdel = cp.random.exponential(TAUAPFAST, nap, dtype='float32') + cp.random.exponential(TAUAPSLOW, nap, dtype='float32')
        # Select wich signals will have ap
        apSig = randint(0, n, dtype='int32')
        tap = (apdel / sampling).astype('int32') + t[apSig]
        hap = 1 - cp.exp(-apdel / TFALL)  # Pulse height as RC circuit
        hap = hap[:, None].astype('float32')
        gainvar = cp.random.normal(1, CCGV, (nap, 1), dtype='float32')
        vect = (cp.arange(SIGPTS, dtype='int32') + cp.zeros((nap, 1), dtype='int32') - tap[:, None])
        vect[vect < 0] = 0
        sig += normpe * cp.sum(signalShapeGPU(vect, TFALL, TRISE, gainvar, hap), axis=0)
    return sig


def PulseCPU(t, h, gainvar, nap):
    """
    PulseCPU(t,h)

    Function that generates the signal from a single SiPM cell.
    This is the "full" version that computes the signal shape on CPU by evaluating the signal shape function.

    Parameters
    ----------
    t : int32
            Time at which the cell is triggered
    h : float32
            The relative pulse height of the cell signal
    gainvar : float32
            Value of cell to cell gain variation for this signal
    nap : int32
            Number of afterpulses in this signal


    Returns
    -------
    s : np.ndarray
            Array containing the generated cell signal
    """

    sig = signalgenfortran(t, h, TFALL, TRISE, SIGPTS, gainvar)    # Calculate signal
    if nap > 0:  # If there are afterpulses generate theyr signals
        for _ in range(nap):
            # APs have a time delay exponential distribution
            apdel = random.expovariate(1 / TAUAPFAST) + random.expovariate(1 / TAUAPSLOW)
            tap = np.int32(apdel / SAMPLING + t)
            hap = 1 - exp(-apdel / TFALL)
            sig += signalgenfortran(tap, hap, TFALL, TRISE, SIGPTS, gainvar)
    return sig


# Function that passes signals times and height to main function for generating signals
def SiPMSignalAction(times, sigH, SNR, BASESPREAD):
    """
    signalGen(times,sigH,SNR,BASESPREAD)

    Function that passes signal height and times to the main function that generates signals. Also adds noise.
    If the number of signals to generate is small uses CPU, else uses GPU to speed up the computation.

    Parameters
    ----------
    times : np.ndarray(int32)
    	Array containing the time at wich SiPM cells are fired, including xt events (sorted)
    sigH : np.ndarray(float32)
    	Array containing the pulse height of each fired SiPM cell
    SNR : double
    	The signal to noise ratio of the noise to add
    BASESPREAD : double
    	Sigma of the value to add as baseline

    Returns
    --------
    signal : np.ndarray
        Array containing the generated SiPM signal
    """

    baseline = random.gauss(0, BASESPREAD)
    if (times.size < CPUTHRESHOLD) or (times.size > GPUMAX):
        signal = np.random.normal(baseline, SNR, SIGPTS)
        gainvars = np.random.normal(1, CCGV, size=times.size)   # Each signal has a ccgv
        naps = poisson(AP, size=times.size)   # Generate number of afterpulses
        for i in range(times.size):
            signal += PulseCPU(times[i], sigH[i], gainvars[i], naps[i])
        return signal
    else:
        signal = cp.random.normal(baseline, SNR,SIGPTS,dtype='float32')
        signal += PulseGPU(cp.asarray(times), cp.asarray(sigH))
        return cp.asnumpy(signal)
