# In this file I define all the functions I will use in the main file of simulation
from libs.FortranFunctions import rollfortran, signalgenfortran, sortfortran
from libs.FortranFunctions import frandom
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
            apdel = frandom.randexp(TAUAPFAST, 1) + frandom.randexp(TAUAPSLOW, 1)
            tap = np.int32(apdel / SAMPLING + t)
            hap = 1 - exp(-apdel / TFALL)
            sig += signalgenfortran(tap, hap, TFALL, TRISE, SIGPTS, gainvar)
    return sig


# Function that passes signals times and height to main function for generating signals
def SiPMSignalAction(times, sigH, SNR, BASESPREAD):
    """
    signalGen(times,sigH,SNR,BASESPREAD)

    Function that passes signal height and times to the main function that
    generates single signals. Also adds noise.

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

    baseline = random.gauss(0, BASESPREAD)  # Add a baseline to the signal
    # Start with gaussian noise
    signal = frandom.randn(baseline, SNR, SIGPTS)
    gainvars = frandom.randn(1, CCGV, times.size)   # Each signal has a ccgv
    naps = frandom.randpoiss(AP, times.size)   # Generate number of afterpulses
    for i in range(times.size):
        signal += PulseCPU(times[i], sigH[i], gainvars[i], naps[i])
    return(signal)
