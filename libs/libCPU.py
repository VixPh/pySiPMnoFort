# In this file I define all the functions I will use in the main file of simulation
from libs.FortranFunctions import signalgenfortran
from variables import *

###############################################################################
##>>>   EDITING THIS FILE MAY SERIOUSLY COMPROMISE SIMULATION BEHAVIOUR   <<<##
###############################################################################
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
    SiPMSignalAction(times,sigH,SNR,basespread)

    Function that passes signal height and times to the main function that generates single signals.
    Also adds noise.

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

    baseline = random.gauss(0, BASESPREAD)  # Add a baseline to the signal
    # Start with gaussian noise
    signal = np.random.normal(baseline, SNR, SIGPTS)
    gainvars = np.random.normal(1, CCGV, size=times.size)   # Each signal has a ccgv
    naps = poisson(AP,size=times.size)   # Generate number of afterpulses
    for i in range(times.size):
        signal += PulseCPU(times[i], sigH[i], gainvars[i], naps[i])
    return(signal)
