# In this file I define all the functions I will use in the main file of simulation
from libs.FortranFunctions import rollfortran, signalgenfortran
from variables import *

###############################################################################
##>>>   EDITING THIS FILE MAY SERIOUSLY COMPROMISE SIMULATION BEHAVIOUR   <<<##
###############################################################################

###GENERATION OF DCR TIMES###
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
        ndcr = poisson(rate * SIGLEN * 1e-9)  # Poissonian distribution
        # Times uniformly distributed
        dcrTime = np.random.random(ndcr) * SIGLEN
    else:
        dcrTime = [0]
        while dcrTime[-1] < SIGLEN:  # Generate from exponential distribution of delays
            delayDCR = random.expovariate(rate) * 1e9
            dcrTime.append(dcrTime[-1] + delayDCR)
        dcrTime = np.array(dcrTime[1:-1])

    return dcrTime


###FUNCTION TO UPDATE THE SIPM MATRIX ARRAY###
def SiPMEventAction(time, XT):
    """
    evtsGen(time,xt)

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
    """
    idx = -1
    evtTimes = []
    sigHtemp = []
    time = np.delete(time, np.where((time > SIGLEN) | (time < 0)))
    n = time.size
    if n > 0:
        idx = randint(0, NCELL, n, dtype=np.int16)
        time = (sort(time) / SAMPLING)
        if FASTXT:
            nxt = poisson(XT * n)
            if nxt > 0:
                addcells = (-1, 1, -CELLSIDE, CELLSIDE, -1 - CELLSIDE,
                            -1 + CELLSIDE, 1 - CELLSIDE, 1 + CELLSIDE)

                xtidx = randint(0, n, nxt, dtype=np.int16)
                xtcells = idx[xtidx]
                xttimes = time[xtidx]
                xtcells += choice(addcells, nxt)
                idx = hstack((idx, xtcells))
                time = hstack((time, xttimes))
        elif not FASTXT:
            addcells = (-1, 1, -CELLSIDE, CELLSIDE, -1 - CELLSIDE,
                        -1 + CELLSIDE, 1 - CELLSIDE, 1 + CELLSIDE)

            nxt = poisson(XT, n)
            if np.sum(nxt) > 0:
                xtcells = []
                xttimes = []
                for i in range(n):
                    if nxt[i] > 0:
                        xtcells.extend(idx[i] + choice(addcells, nxt[i]))
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
                mCellT = sort(time[idx == mCellI])
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

    if not args.Graphics:
        idx = None
    return(np.array(evtTimes, dtype='int32'), np.array(sigH, dtype='float32'), idx)


### GENERATION OF SIGNALS SHAPES FAST ###
if args.signal is None:  # If generating signals fast (default)
    x = np.arange(0, SIGPTS)
# Define the model of my signal (calculate it only once)
    signalmodel = signalgenfortran(x, 1, TFALL, TRISE, SIGPTS, 1)
    def SiPMSignalAction(times, sigH, SNR, BASESPREAD):
        """
        signalGen(times,sigH,SNR,basespread)

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
        basespread : double
                Sigma of the value to add as baseline

        Returns
        --------
        signal : np.ndarray
                Array containing the generated SiPM signal
        """
        baseline = random.gauss(0, BASESPREAD)  # Add a baseline
        signal = np.random.normal(baseline, SNR, SIGPTS)    # Start with gaussian noise
        for i in range(times.size):   # Generate signals and sum them
            signal += PulseCPU(times[i], sigH[i])
        return signal

    def PulseCPU(t, h):
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

        Returns
        -------
        s : np.ndarray
                Array containing the generated cell signal
        """
        gainvar = random.gauss(1, CCGV)    # Each signal has a ccgv
        sig = rollfortran(signalmodel, t, gainvar, h, SIGPTS)   # Move the model signal
        nap = poisson(AP)  # Generate number of afterpulses
        if nap > 0:
            for _ in range(nap):
                # APs have a time delay exponential distribution
                apdel = random.expovariate(1 / TAUAPFAST) + random.expovariate(1 / TAUAPSLOW)
                tap = np.int32(apdel / SAMPLING + t)
                if tap < SIGPTS:
                    # Generate ap signal height as an RC circuit
                    hap = 1 - exp(-apdel / CELLRECOVERY)
                    sap = rollfortran(signalmodel, tap, gainvar, hap, SIGPTS)
                    sig += sap       # Add each ap
        return sig


### FULL GENERATION OF SIGNAL SHAPES ###
else:								# Recalculate the signal each time
    if args.device == 'cpu':
        from libs.libCPU import *   # Generation fully on cpu
    if args.device == 'gpu':
        from libs.libGPU import *   # Cpu for low light and cpu for high light

###SOME STATISCTICS AT END OF SCRIPT###


def somestats(output):
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
    other = output[:, 5]

    ROOT.gROOT.SetStyle("ATLAS")
    ROOT.gStyle.SetOptStat(1)
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetHistFillColor(603)
    c1 = ROOT.TCanvas('c1', 'Histogram', 600, 400)
    c2 = ROOT.TCanvas('c2', 'Histogram', 600, 400)
    c3 = ROOT.TCanvas('c3', 'Histogram', 600, 400)
    c4 = ROOT.TCanvas('c4', 'Histogram', 600, 400)
    c5 = ROOT.TCanvas('c5', 'Staircase', 600, 400)
    c6 = ROOT.TCanvas('c6', 'Histogram', 600, 400)
    c7 = ROOT.TCanvas('c7', 'Histogram', 600, 400)

    h1 = ROOT.TH1F('Integral', 'Integral', 750, integral.min(), integral.max())
    h1.SetXTitle('Integral [A.U.]')
    h1.SetYTitle('Entries')
    h2 = ROOT.TH1F('Peak Value', 'Peak', 750, peak.min(), peak.max())
    h2.SetXTitle('Peak [A.U.]')
    h2.SetYTitle('Entries')
    inf = tstart.min()
    sup = tstart.max()
    n_bin = int((sup - inf) / SAMPLING)
    h3 = ROOT.TH1F('Starting Time', 'Tstart', n_bin, inf, sup)
    h3.SetXTitle('Starting Time [ns]')
    h4 = ROOT.TH2F('Histogram', 'Peak - Integral', 200, min(peak),max(peak), 200, min(integral), max(integral))
    h4.SetXTitle('Peak [A.U.]')
    h4.SetYTitle('Integral [A.U.]')
    inf = 0
    sup = tovert.max()
    n_bin = int((sup - inf) / SAMPLING)
    h5 = ROOT.TH1F('ToT', 'Time over threshld', n_bin, inf, sup)
    h5.SetXTitle('Time [ns]')
    h5.SetYTitle('Entries')
    inf = 0
    sup = tpeak.max()
    n_bin = int((sup - inf) / SAMPLING)
    h6 = ROOT.TH1F('ToP', 'Time of peak', n_bin, inf, sup)
    h6.SetXTitle('Time [ns]')
    h6.SetYTitle('Entries')
    [h1.Fill(i) for i in integral]
    [h2.Fill(i) for i in peak]
    [h3.Fill(i) for i in tstart]
    [h4.Fill(i, j) for i, j in zip(peak, integral)]
    [h5.Fill(i) for i in tovert]
    [h6.Fill(i) for i in tpeak]
    c1.cd()
    h1.Draw("bar9")
    c2.cd()
    h2.Draw("bar9")
    c3.cd()
    h3.Draw("E9")
    c4.cd()
    h4.Draw("colz")
    c5.cd()
    c5.SetLogy()
    staircase = h2.GetCumulative(False, ' Staircase')
    staircase.Scale(1e-3 / h2.Integral() / (INTGATE * SAMPLING * 1e-9))
    staircase.SetLineColor(2)
    staircase.SetTitle('Staircase')
    staircase.SetXTitle('Peak [A.U.]')
    staircase.SetYTitle('Rate [kHz]')
    staircase.SetFillColor(0)
    staircase.Draw("hist9")
    c6.cd()
    h5.Draw("bar9")
    c7.cd()
    h6.Draw("bar9")
    c1.Update()
    c2.Update()
    c3.Update()
    c4.Update()
    c5.Update()
    c6.Update()
    c7.Update()
    input('Press <RET> to continue...')
    return


drawn = [False] * nJobs
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

    sipmmatrix = np.zeros((CELLSIDE, CELLSIDE), dtype='float32')

    if np.all(idx >= 0):
        times = np.hstack((sigTimes, dcrTime))
        rows = idx // CELLSIDE
        cols = idx % CELLSIDE
        sipmmatrix[rows, cols] = times

    if not drawn[current_core]:
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
        plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99,
                            top=0.99, wspace=0.05, hspace=0.15)
    else:
        axlist = plt.gcf().axes
        ax = axlist[0]
        axm = axlist[1]

    textstring = f"Core: {current_core:d}\n"
    textstring += f"Device: {dev:s}\n"
    textstring += f"Photons:{sigTimes.size-dcrTime.size:d} Dcr:{dcrTime.size:d}\n"
    txt = ax.text(0.5, 0.70, textstring, transform=ax.transAxes, fontsize=10)

    ax.plot(np.arange(SIGPTS) * SAMPLING, signal, '-b', linewidth=0.5)
    axm.matshow(sipmmatrix)

    plt.pause(0.25)
    ax.lines[-1].remove()
    txt.remove()


def initializeRandomPool():
    """
    initializeRandomPool()

    Function that initializes random seeds for each worker in the multiprocessing Pool
    """
    # Get current core number (On MacOs psutil is not working)
    current_core = multiprocessing.current_process().name
    current_core = int(current_core.split('-')[-1])
    time.sleep(0.5 / nJobs * current_core)
    # Get some random bits from the sistem entropy pool
    rngseed = struct.unpack('I', os.urandom(4))[0] + current_core
    random.seed(rngseed)   # Change rng seed for each worker
    np.random.seed(rngseed)
    core = multiprocessing.current_process().name
    print("Initializing simulation on %s with seed %d\r" % (core, rngseed))


output = []
def Callback(results):
    # Function that saves results from multiprocessing pool
    output.append(results)


def SaveFile(fname, out):
    integral = out[:, 0]
    peak = out[:, 1]
    tstart = out[:, 2]
    tover = out[:, 3]
    ptime = out[:, 4]
    other = np.vstack(out[:, 5])

    event = other[:, 0]
    fiber = other[:, 1]
    theta = other[:, 2]
    phi = other[:, 3]

    f = uproot.recreate(fname)
    f['SiPMData'] = uproot.newtree(
        {'Integral': 'float32',
         'Peak': 'float32',
         'ToA': 'float32',
         'ToT': 'float32',
         'ToP': 'float32'})
    f['GeometryData'] = uproot.newtree(
        {'EventId': 'int32',
         'FiberId': 'int32',
         'FiberTheta': 'float32',
         'FiberPhi': 'float32'})

    f['SiPMData']['Integral'].newbasket(integral)
    f['SiPMData']['Peak'].newbasket(peak)
    f['SiPMData']['ToA'].newbasket(tstart)
    f['SiPMData']['ToP'].newbasket(ptime)
    f['SiPMData']['ToT'].newbasket(tover)

    f['GeometryData']['EventId'].newbasket(event)
    f['GeometryData']['FiberId'].newbasket(fiber)
    f['GeometryData']['FiberTheta'].newbasket(theta)
    f['GeometryData']['FiberPhi'].newbasket(phi)


def SaveWaves(fname, out):
    signals = np.vstack(out[:, -1])

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

    with h5py.File('waveforms.hdf5', 'a') as hf:
        fname = fname + '0'
        groups = list(hf.keys())
        if fname in groups:
            i = int(groups[-1][-1])
            fname = fname[:-1] + str(i + 1)

        grp = hf.create_group(fname)
        dset1 = grp.create_dataset('Waveforms',
                                   shape=(signals.shape),
                                   dtype='f',
                                   compression='lzf',
                                   chunks=(1, signals.shape[1]))
        dset2 = grp.create_dataset('SiPMSettings',
                                   shape=(len(sipmsettings),),
                                   dtype='f',
                                   compression='gzip',
                                   compression_opts=9)
        dset1[...] = signals
        dset2[...] = sipmsettings
