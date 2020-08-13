# In this file I define all the functions I will use in the main file of simulation
from libs.FortranFunctions import rollfortran, signalgenfortran, sortfortran
from libs.FortranFunctions import frandom
from variables import *


###############################################################################
##>>>   EDITING THIS FILE MAY SERIOUSLY COMPROMISE SIMULATION BEHAVIOUR   <<<##
###############################################################################

### GENERATION OF DCR TIMES ###
def addDCR(rate):
    """
    addDCR(rate)

    Function that generates times of dcr events.

    Parameters
    ------------
    rate : double
            Rate of dcr in kHz

    Returns
    --------
    dcrTime : np.ndarray
            Array containing dcr times
    """

    if FASTDCR:
        ndcr = frandom.randpoiss(rate * SIGLEN * 1e-9, 1)  # Poissonian distribution
        # Times uniformly distributed
        dcrTime = np.random.random(ndcr) * SIGLEN
    else:
        dcrTime = [0]
        while dcrTime[-1] < SIGLEN:  # Generate from exponential distribution of delays
            delayDCR = frandom.randexp(1/rate, 1) * 1e9
            dcrTime.extend(dcrTime[-1] + delayDCR)
        dcrTime = np.array(dcrTime[1:-1])

    return dcrTime


### FUNCTION TO UPDATE THE SIPM MATRIX ARRAY ###
def SiPMEventAction(time, XT):
    """
    evtsGen(time,XT)

    Function that calculates signal height for each SiPM cell and adds optical crosstalk events.

    Parameters
    ----------
    time : np.ndarray
        Array containing the time at which SiPM cells are fired
    XT : double
        Value of optical crosstalk probability

    Returns
    -------
    evtTimes : np.ndarray(int32)
            Array containing the time at wich SiPM cells are fired, including xt events (sorted)
    sigH : np.ndarray(float32)
            Array containing the pulse height of each fired SiPM cell
    idx : np.ndarray(int16)
            Array containing the unique index of the hitted cells
    """

    idx = -1
    evtTimes = []
    sigHtemp = []
    if np.any(time > SIGLEN):
        print('Detected events past the signal length, deleting them...')
        time = time[time < SIGLEN]
    n = time.size
    if n > 0:
        idx = frandom.randint(NCELL, n)
        sortfortran(time)
        time /= SAMPLING
        addcells = np.array((-1, 1, -CELLSIDE, CELLSIDE, -1 - CELLSIDE, -1 + CELLSIDE, 1 - CELLSIDE, 1 + CELLSIDE))
        if FASTXT:
            nxt = frandom.randpoiss(XT * n, 1)
            if nxt:
                xtidx = frandom.randint(n-1, nxt)
                xtcells = idx[xtidx]
                xttimes = time[xtidx]
                xtcells += addcells[frandom.randint(7, nxt)]
                idx = hstack((idx, xtcells))
                time = hstack((time, xttimes))
        else:
            nxt = frandom.randpoiss(XT, n)
            if np.count_nonzero(nxt):
                xtcells = []
                xttimes = []
                for i in range(n):
                    if nxt[i]:
                        xtcells.extend(idx[i] + addcells[frandom.randint(7, nxt[i])])
                        xttimes.extend([time[i]] * nxt[i])
                idx = hstack((idx, xtcells))
                time = hstack((time, xttimes))

        if idx.size == len(set(idx)):
            evtTimes = time
            sigH = time * 0 + 1
        else:
            _, uCellsIdx, uCellsCounts = np.unique(idx, return_index=True, return_counts=True)
            # Times of cells fired once
            sTimes = time[uCellsIdx[uCellsCounts == 1]]
            evtTimes.append(sTimes)
            sigHtemp.append(sTimes * 0 + 1)

            # Idx of cell fired multple times
            mCellIdx = idx[uCellsIdx[uCellsCounts > 1]]
            for i in range(mCellIdx.size):
                mCellI = mCellIdx[i]    # Loop on cells fired multiple times
                # Times of events in the same cell
                mCellT = time[idx == mCellI]
                sortfortran(mCellT)
                # Delays of events consecutive to the first
                delays = mCellT - mCellT[0]
                # Calculate height of pulses
                h = 1 - exp(-delays / CELLRECOVERY)
                h[0] = 1
                evtTimes.append(mCellT)
                sigHtemp.append(h)
            evtTimes = hstack(evtTimes)
            sigH = hstack(sigHtemp)
    else:
        # If signal has 0 pe and 0 dcr pass an empty array
        evtTimes = np.empty(0)
        sigH = np.empty(0)

    return(np.int32(evtTimes), np.float32(sigH), idx)


### GENERATION OF SIGNALS SHAPES FAST ###
if args.signal is None:  # If generating signals fast (default)
    x = np.arange(0, SIGPTS)
# Define the model of my signal (calculate it only once)
    signalmodel = signalgenfortran(0, 1, TFALL, TRISE, SIGPTS, 1)

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
        baseline = random.gauss(0, BASESPREAD)  # Add a baseline
        signal = frandom.randn(baseline, SNR, SIGPTS)    # Start with gaussian noise
        gainvars = frandom.randn(1, CCGV, times.size)    # Each signal has a ccgv
        naps = frandom.randpoiss(AP, times.size)
        for i in range(times.size):   # Generate signals and sum them
            signal += PulseCPU(times[i], sigH[i], gainvars[i], naps[i])
        return signal

    def PulseCPU(t, h, gainvar, nap):
        """
        PulseCPU(t,h)

        Function that generates the signal from a single SiPM cell.
        This is the "fast" version that uses a pre-computed signal
        shape and translates it in time.

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

        sig = rollfortran(signalmodel, t, gainvar, h)   # Move the model signal
        if nap:
            for _ in range(nap):
                # APs have a time delay exponential distribution
                apdel = frandom.randexp(TAUAPFAST, 1) + frandom.randexp(TAUAPSLOW, 1)
                tap = np.int32(apdel / SAMPLING + t)
                if tap < SIGPTS:
                    # Generate ap signal height as an RC circuit
                    hap = 1 - exp(-apdel / CELLRECOVERY)
                    sap = rollfortran(signalmodel, tap, gainvar, hap)
                    sig += sap       # Add each ap
        return sig


### FULL GENERATION OF SIGNAL SHAPES ###
else:								# Recalculate the signal each time
    if args.device == 'cpu':
        from libs.libCPU import *   # Generation fully on cpu
    if args.device == 'gpu':
        from libs.libGPU import *   # Cpu for low light and cpu for high light

### SOME STATISCTICS AT END OF SCRIPT ###


def somestats(output, realpe=None):
    """
    somestats(output)

    Function that displays histograms of generated events

    Parameters
    ----------
    output : np.ndarray
            Array containing output of the simulation.
            This array contains the integral,
            peak and starting time in the first three columns.

    Note
    -----
    See docs for a detailed description of integral, peak and tstart
    """

    integral = output[:, 0]
    peak = output[:, 1]
    tstart = output[:, 2]
    tovert = output[:, 3]
    tpeak = output[:, 4]

    plt.figure()
    plt.title('Integral')
    plt.hist(integral, 500, color='k')
    plt.xlabel('Integrated charge [A.U.]')

    plt.figure()
    plt.title('Peak')
    plt.hist(peak, 500, color='k')
    plt.xlabel('Peak value [A.U.]')

    plt.figure()
    plt.title('Starting time')
    plt.hist(tstart, np.arange(tstart.min()-2, tstart.max()+2, 2*SAMPLING), color='k')
    plt.xlabel('Starting time [ns]')
    plt.yscale('log')

    plt.figure()
    plt.title('Time over threshold')
    plt.hist(tovert, np.arange(tovert.min(), tovert.max(), 2*SAMPLING), color='k')
    plt.xlabel('Time over threshold [ns]')

    plt.figure()
    plt.title('Peaking time')
    plt.hist(tpeak, np.arange(tpeak.min(), tpeak.max(), 2*SAMPLING), color='k')
    plt.xlabel('Peaking time [ns]')

    plt.figure()
    x = np.sort(peak)
    y = np.empty_like(x)
    for i in range(x.size):
        y[i] = np.count_nonzero(x > x[i])

    y *= (1e-3 / y.size / (INTGATE * SAMPLING * 1e-9))
    if peak.max() < 5:
        plt.hlines(DCR/1e3, x.min(), x.max(), 'k', label=f'DCR = {DCR*1e-3:.0f} kHz')
        plt.legend()
    plt.plot(x, y, '.r')
    plt.yscale('log')
    plt.ylabel('Rate [kHz]')
    plt.xlabel('Threshold')
    plt.title('Staircase')

    if realpe is not None:
        plt.figure()
        plt.subplot(121)
        plt.scatter(realpe, peak, c='k', s=2)
        plt.xlabel('Real number of photoelectrons')
        plt.ylabel('Peak value measured')
        plt.subplot(122)
        plt.hist2d(realpe, peak, bins=(100, 100))
        plt.xlabel('Real number of photoelectrons')
        plt.ylabel('Peak value measured')

    plt.figure()
    plt.subplot(121)
    plt.scatter(peak, integral, c='k', s=2)
    plt.xlabel('Peak value measured')
    plt.ylabel('Integrated charge')
    plt.subplot(122)
    plt.hist2d(peak, integral, bins=(100, 100))
    plt.xlabel('Peak value measured')
    plt.ylabel('Integrated charge')
    plt.show()


drawn = [False] * nJobs
opened = [True] * nJobs
def sigPlot(signal, sigTimes, dcrTime, dev, idx):
    """
    sigPlot(signal,sigTimes,dcrTime,dev)

    Function that plots each signal pulse produced in the simulation.

    Parameters
    -----------
    signal : np.ndarray
            Array containing the generated SiPM signal
    sigTimes : np.ndarray
            Array containing all photons events times (including xt)
    dcrTime : np.ndarray
            Array containing dcr events times
    dev : str
            String that describes the device on which the signal is computed
    """

    current_core = multiprocessing.current_process().name
    if current_core == 'MainProcess':
        current_core = 0
    else:
        current_core = int(current_core.split('-')[-1]) - 1

    sipmmatrix = np.zeros((CELLSIDE, CELLSIDE), dtype='bool')

    if np.any(idx >= 0):
        times = sigTimes / SAMPLING
        rows = (idx // CELLSIDE) % CELLSIDE
        cols = (idx % CELLSIDE) % CELLSIDE
        sipmmatrix[rows, cols] = True

    textstring = f"Core: {current_core:d}\n"
    textstring += f"Device: {dev:s}\n"
    textstring += f"Photons:{sigTimes.size-dcrTime.size:d} Dcr:{dcrTime.size:d}\n"
    if not drawn[current_core]:
        timearray = np.arange(SIGPTS) * SAMPLING
        ax = plt.subplot(211)
        screenx = 1920
        screeny = 1080
        xsize = int(screenx / 6)
        ysize = int(screeny / 2)
        xpos = (current_core % 6) * xsize
        ypos = ((current_core // 6) % 3) * ysize
        plt.get_current_fig_manager().window.setGeometry(xpos, ypos, xsize, ysize)
        ax.hlines(-0.5, 0, INTSTART * SAMPLING, 'r')
        ax.vlines(INTSTART * SAMPLING, -0.5, -1, 'r')
        ax.vlines((INTSTART - PREG) * SAMPLING, -0.5, -1, 'r')
        ax.hlines(-1, INTSTART * SAMPLING, (INTSTART + INTGATE) * SAMPLING, 'r')
        ax.vlines((INTSTART + INTGATE) * SAMPLING, -1, -0.5, 'r')
        ax.hlines(-0.5, (INTSTART + INTGATE) * SAMPLING, SIGLEN, 'r')
        ax.grid(linestyle=':')
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        drawn[current_core] = True

        axm = plt.subplot(212)
        axm.set_axis_off()
        plt.subplots_adjust(left=0.1, bottom=0.01, right=0.99,
                            top=0.99, wspace=0.05, hspace=0.15)
        ax.plot(timearray, signal, '-b', linewidth=0.5)
        axm.matshow(sipmmatrix, aspect='equal', filternorm=False, resample=False, cmap='binary_r')
    else:
        axlist = plt.gcf().axes
        if len(axlist) == 0:
            opened[current_core] = False
        if opened[current_core]:
            ax = axlist[0]
            axm = axlist[1]
            if signal.max() > opened[current_core]:
                opened[current_core] = True

    if opened[current_core]:
        line = ax.lines[-1]
        img = axm.get_images()[0]

        txt = ax.text(0.5, 0.70, textstring, transform=ax.transAxes, fontsize=10)

        line.set_ydata(signal)
        ax.relim()
        ax.autoscale_view(tight=False, scalex=False, scaley=True)
        img.set_data(sipmmatrix)

        plt.pause(float(args.Graphics)/1000)
        txt.remove()


def initializeRandomPool():
    """
    initializeRandomPool()

    Function that initializes random seeds for each worker in the multiprocessing Pool
    """

    # Get current core number (On MacOs psutil is not working)
    current_core = multiprocessing.current_process().name
    current_core = int(current_core.split('-')[-1])
    # Get some random bits from the sistem entropy pool
    rngseed = struct.unpack('I', os.urandom(4))[0] + current_core
    random.seed(rngseed)   # Change rng seed for each worker
    np.random.seed(rngseed)
    core = multiprocessing.current_process().name
    print("Initializing simulation on %s with seed %d\r" % (core, rngseed))


def SaveFile(fname, out, other=None):
    integral = out[:, 0]
    peak = out[:, 1]
    tstart = out[:, 2]
    tover = out[:, 3]
    ptime = out[:, 4]

    f = uproot.recreate(fname)
    f['SiPMData'] = uproot.newtree(
        {'Integral': 'float32',
         'Peak': 'float32',
         'ToA': 'float32',
         'ToT': 'float32',
         'ToP': 'float32'})

    f['SiPMData']['Integral'].newbasket(integral)
    f['SiPMData']['Peak'].newbasket(peak)
    f['SiPMData']['ToA'].newbasket(tstart)
    f['SiPMData']['ToP'].newbasket(ptime)
    f['SiPMData']['ToT'].newbasket(tover)

    if other:
        other = np.array(other)
        event = np.int32(other[:, 0])
        fibertype = other[:, 1] == 'Scin'
        fiberid = np.int64(other[:, 2])
        x = np.float32(other[:, 3])
        y = np.float32(other[:, 4])
        z = np.float32(other[:, 5])

        f['GeometryData'] = uproot.newtree(
            {'EventId': np.int32,
             'FiberType': np.int8,
             'FiberId': np.int64,
             'FiberX': np.float32,
             'FiberY': np.float32,
             'FiberZ': np.float32})

        f['GeometryData']['EventId'].newbasket(event)
        f['GeometryData']['FiberType'].newbasket(fibertype)
        f['GeometryData']['FiberId'].newbasket(fiberid)
        f['GeometryData']['FiberX'].newbasket(x)
        f['GeometryData']['FiberY'].newbasket(y)
        f['GeometryData']['FiberZ'].newbasket(z)



def SaveWaves(fname, signals):
    while os.path.exists(fname):
        f = fname.split('.')
        if f[0][-1].isnumeric():
            fname = f[0][:-1]+str(int(f[0][-1])+1)+'.hdf5'
        else:
            fname = f[0]+'0'+'.'+f[1]

    sipmsettings = [SIZE,
                    CELLSIZE,
                    DCR,
                    XT,
                    AP,
                    CELLRECOVERY * SAMPLING,
                    TRISE * SAMPLING,
                    TFALL * SAMPLING,
                    -20 * np.log10(SNR**2),
                    CCGV]

    with h5py.File(fname, 'a') as hf:
        dset1 = hf.create_dataset('Waveforms',
                                  shape=(signals.shape),
                                  dtype='f',
                                  compression='gzip',
                                  chunks=(1, signals.shape[1]),
                                  compression_opts=9)
        dset2 = hf.create_dataset('SiPMSettings',
                                  shape=(len(sipmsettings),),
                                  dtype='f',
                                  compression='gzip',
                                  compression_opts=9)
        dset1[...] = signals
        dset2[...] = sipmsettings
