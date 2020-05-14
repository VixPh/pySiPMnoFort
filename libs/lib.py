#In this file I define all the functions I will use in the main file of simulation
from variables import *
from libs.FortranFunctions import signalgenfortran
from libs.FortranFunctions import rollfortran
#################################################################################################################
################>>>   EDITING THIS FILE MAY SERIOUSLY COMPROMISE SIMULATION BEHAVIOUR   <<<######################
#################################################################################################################

###GENERATION OF DCR TIMES AND INDICES###
def addDCR(rate):											# Generate dcr times
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
	if fastDCR:
		ndcr = poisson(rate*siglen*1e-9)					# Number of dcr events folows a poissonian distribution
		dcrTime = np.random.random(ndcr)*siglen			# dcr times are uniformly distributed
	else:
		dcrTime = [0]
		while dcrTime[-1]<siglen:						# Generate dcr starting from exponential distribution of delays
			delayDCR = random.expovariate(rate)*1e9
			dcrTime.append(dcrTime[-1]+delayDCR)
		dcrTime = np.array(dcrTime[1:-1])
	return(dcrTime)

###FUNCTION TO UPDATE THE SIPM MATRIX ARRAY###
def SiPMEventAction(time,xt):
	"""
	evtsGen(time,xt)

	Function that calculates signal height for each SiPM cell and adds optical crosstalk events.

	Parameters
	----------
	time : np.ndarray
		Array containing the time at which SiPM cells are fired
	xt : double
		Value of optical crosstalk probability

	Returns
	-------
	evtTimes : np.ndarray(int32)
		Array containing the time at wich SiPM cells are fired, including xt events (sorted)
	sigH : np.ndarray(float32)
		Array containing the pulse height of each fired SiPM cell
	"""
	evtTimes = []
	sigHtemp = []
	time = np.delete(time,np.where((time>siglen)|(time<0)))
	n = time.size															# Number of hitted cells
	if n > 0:
		idx = randint(0,ncell,n,dtype=np.int16)								# Generate random indexes for each cell to be fired (suppose uniform distribution on sensor surface)
		time = (sort(time)/sampling)										# Convert times in sampling times units and sort them
		if fastXT:
			nxt = poisson(xt*n)												# Generate number of cross talk events (xt per each cell so n*xt for n cells)
			if nxt > 0:
				addcells=(-1,1,-cellside,cellside,-1-cellside,-1+cellside,1-cellside,1+cellside)	#Values to add to change cell for XT (8 adjacent cells)
				xtidx = randint(0,n,nxt,dtype=np.int16)						# Indexes of cells triggering xt events
				xtcells = idx[xtidx]										# Cell ids triggering xt events
				xttimes = time[xtidx]										# Times of xt events (same of the cell that generates it)
				xtcells += choice(addcells,nxt)								# Select neighbouring cells
				idx = hstack((idx,xtcells))									# Add XT cells to the list
				time = hstack((time,xttimes))								# Add XT times to the list
		elif not fastXT:
			addcells=(-1,1,-cellside,cellside,-1-cellside,-1+cellside,1-cellside,1+cellside)
			nxt = poisson(xt,n)												# Generate number of poisson events per each cell
			if np.sum(nxt)>0:												# If I have any
				xtcells = []
				xttimes = []
				for i in range(n):											# Per each cell add xt as above
					if nxt[i]>0:
						xtcells.extend(idx[i]+choice(addcells,nxt[i]))
						xttimes.extend([time[i]]*nxt[i])
				idx = hstack((idx,xtcells))
				time = hstack((time,xttimes))

		if idx.size == len(set(idx)):											# If all cells index are different output all times without other iterations
			evtTimes = time
			sigH = time*0+1														# All signals have height 1
		else:
			uCells,uCellsIdx,uCellsCounts = np.unique(idx,return_index=True,return_counts=True)
			sTimes = time[uCellsIdx[uCellsCounts==1]]							# Times of cells fired once (nothing to do)
			evtTimes.append(sTimes)												# Add times of cells fired once
			sigHtemp.append(sTimes*0+1)											# Cells fired once have an height of 1

			mCellIdx = idx[uCellsIdx[uCellsCounts>1]]							# Idx of cell fired multple times (consider cell recovery)
			for i in range(mCellIdx.size):
				mCellI = mCellIdx[i]											# Loop on cells fired multiple times
				mCellT = sort(time[idx==mCellI])								# Times of events in the same cell # TODO: Do I need to sort again?
				delays = mCellT-mCellT[0]										# Delays of events consecutive to the first one (first delay returns 0)
				h = 1-exp(-delays/cellrecovery)									# Calculate height of pulses (RC discharge circuit)
				h[0] = 1														# First pulse height is one
				evtTimes.append(mCellT)
				sigHtemp.append(h)
			evtTimes = hstack(evtTimes)											# Stack al times and heights of signals
			sigH = hstack(sigHtemp)
	else:
		evtTimes = np.empty(0)													# If signal has 0 pe and 0 dcr pass an empty array
		sigH = np.empty(0)
	return(np.array(evtTimes,dtype='int32'),np.array(sigH,dtype='float32'))

### GENERATION OF SIGNALS SHAPES FAST ###
if args.signal is None:											# If generating signals fast (default)
	x = np.arange(0,sigpts)
	signalmodel = signalgenfortran(x,1,tfall,trise,sigpts,1)	# Define the model of my signal (calculate it only once)

	def SiPMSignalAction(times,sigH,SNR,basespread):							# Function that passes signals times and height to main function for generating signals
		"""
		signalGen(times,sigH,SNR,basespread)

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
		for i in range(times.size):								# Generate signals and sum them
			signal += PulseCPU(times[i],sigH[i])				# Generate signals
		return(signal)

	def PulseCPU(t,h):
		"""
		PulseCPU(t,h)

		Function that generates the signal from a single SiPM cell. This is the "fast" version that uses a pre-computed signal shape and translates it in time.

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
		gainvar = random.gauss(1,ccgv)							# Each signal has a cell to cell gain variation
		s = rollfortran(signalmodel,t,gainvar,h,sigpts)			# Move the model signal at the appropriate time
		nap = poisson(ap)										# Generate number of afterpulses
		if nap>0:
			for i in range(nap):
				apdel = random.expovariate(1/tauapfast)+random.expovariate(1/tauapslow)		# APs have a time delay exponential distribution
				tap = np.int32(apdel/sampling+t)											# Each afterpulse has a delay from its "main" signal that follows a exp distribution
				if tap<siglen:
					hap = 1-exp(-apdel/tfall)													# Generate ap signal height as an RC circuit
					sap = rollfortran(signalmodel,tap,gainvar,hap,sigpts)									# Generate ap signal as above
					s += sap																	# Add each ap
		return(s)

### FULL GENERATION OF SIGNAL SHAPES ###
else:												# In this case recalculate the signal each time
	if args.device == 'cpu':
		from libs.libCPU import *					# Generation fully on cpu
	if args.device == 'gpu':
		from libs.libGPU import *					# Generation on cpu for low light and cpu for high light

###SOME STATISCTICS AT END OF SCRIPT###
def somestats(output):
	"""
	somestats(output)

	Function that displays histograms of generated events

	Parameters
	----------
	output : np.ndarray
		Array containing output of the simulation. This array contains the integral, peak and starting time in the first three columns.

	Note
	-----
	See docs for a detailed description of integral, peak and tstart
	"""

	integral = output[:,0]
	peak = output[:,1]
	tstart = output[:,2]
	tovert = output[:,3]
	tpeak = output[:,4]
	other = output[:,5]

	ROOT.gROOT.SetStyle("ATLAS")
	ROOT.gStyle.SetOptStat(1)
	ROOT.gStyle.SetOptFit(1)
	ROOT.gStyle.SetHistFillColor(603)
	c1 = ROOT.TCanvas('c1','Histogram',600,400)
	c2 = ROOT.TCanvas('c2','Histogram',600,400)
	c3 = ROOT.TCanvas('c3','Histogram',600,400)
	c4 = ROOT.TCanvas('c4','Histogram',600,400)
	c5 = ROOT.TCanvas('c5','Cumulatives',600,400)
	c6 = ROOT.TCanvas('c6','Histogram',600,400)
	c7 = ROOT.TCanvas('c7','Histogram',600,400)



	h1 = ROOT.TH1F('Integral','Integral',750,integral.min(),integral.max())
	h1.SetXTitle('Integral [A.U.]')
	h1.SetYTitle('Entries')
	h2 = ROOT.TH1F('Peak Value','Peak',750,peak.min(),peak.max())
	h2.SetXTitle('Peak [A.U.]')
	h2.SetYTitle('Entries')
	h3 = ROOT.TH1F('Starting Time','Tstart',intgate,tstart.min()-1,tstart.max())
	h3.SetXTitle('Starting Time [ns]')
	h4 = ROOT.TH2F('Histogram','Peak - Integral',200,min(peak),max(peak),200,min(integral),max(integral))
	h4.SetXTitle('Peak [A.U.]')
	h4.SetYTitle('Integral [A.U.]')
	h5 = ROOT.TH1F('ToT','Time over threshld',intgate,0,tovert.max())
	h5.SetXTitle('Time [ns]')
	h5.SetYTitle('Entries')
	h6 = ROOT.TH1F('ToP','Time of peak',intgate,0,tpeak.max())
	h6.SetXTitle('Time [ns]')
	h6.SetYTitle('Entries')
	[h1.Fill(i) for i in integral]
	[h2.Fill(i) for i in peak]
	[h3.Fill(i) for i in tstart]
	[h4.Fill(i,j) for i,j in zip(peak,integral)]
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
	staircase = h2.GetCumulative(False,' Cumulative')
	staircase.Scale(1e-3/h2.Integral()/(intgate*sampling*1e-9))
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
	c1.Update();c2.Update();c3.Update();c4.Update();c5.Update();c6.Update();c7.Update()
	input('Press <RET> to continue...')
	return

drawn = [False]*nJobs
def sigPlot(signal,sigTimes,dcrTime,dev):
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
	if current_core=='MainProcess':
		current_core = 0
	else:
		current_core = int(current_core.split('-')[-1])-1
	ax = plt.gca()
	if not drawn[current_core]:
		screenx = 1920; screeny = 1080;
		xsize = int(screenx/6); ysize = int(screeny/3);
		xpos = (current_core%6)*xsize
		ypos = ((current_core//6)%3)*ysize
		plt.get_current_fig_manager().window.setGeometry(xpos,ypos,xsize,ysize)
		ax.hlines(-0.5,0,intstart*sampling,'r')
		ax.vlines(intstart*sampling,-0.5,-1,'r')
		ax.vlines((intstart-preg)*sampling,-0.5,-1,'r')
		ax.hlines(-1,intstart*sampling,(intstart+intgate)*sampling,'r')
		ax.vlines((intstart+intgate)*sampling,-1,-0.5,'r')
		ax.hlines(-0.5,(intstart+intgate)*sampling,siglen,'r')
		ax.grid(linestyle=':')
		ax.yaxis.set_major_locator(MaxNLocator(integer=True))
		drawn[current_core] = True
		ax.get_figure().subplots_adjust(left=0.075,right=0.99,top=0.99,bottom=0.075)

	textstring = f"Core: {current_core:d}\n"
	textstring += f"Device: {dev:s}\n"
	textstring += f"Photons:{sigTimes.size-dcrTime.size:d} Dcr:{dcrTime.size:d}\n"
	t = ax.text(0.05,0.70,textstring,transform=ax.transAxes,fontsize=12)

	ax.plot(np.arange(sigpts)*sampling,signal,'-b',linewidth=0.5)

	plt.pause(0.05)
	ax.lines[-1].remove()
	t.remove()

def initializeRandomPool():
	"""
	initializeRandomPool()

	Function that initializes random seeds for each worker in the multiprocessing Pool
	"""
	current_core = multiprocessing.current_process().name					# Get current core number (using this trick becouse on MacOs psutil is not working)
	current_core = int(current_core.split('-')[-1])
	time.sleep(0.5/nJobs*current_core)
	rngseed=struct.unpack('I',os.urandom(4))[0]+current_core				# Get some random bits from the sistem entropy pool
	random.seed(rngseed)													# Change rng seed (must be different for each worker, otherwise they will generate same values)
	np.random.seed(rngseed)
	core=multiprocessing.current_process().name
	print("Initializing simulation on %s with seed %d\r"%(core,rngseed))

output = []
def Callback(results):
	output.append(results)													# Function that saves results from multiprocessing pool

def SaveFile(fname,output):
	integral = output[:,0]
	peak = output[:,1]
	tstart = output[:,2]
	tover = output[:,3]
	ptime = output[:,4]
	other = np.vstack(output[:,5])

	event = other[:,0]
	fiber = other[:,1]
	theta = other[:,2]
	phi = other[:,3]

	f = uproot.recreate(fname)
	f['SiPMData'] = uproot.newtree({'Integral':'float32','Peak':'float32','ToA':'float32','ToT':'float32','ToP':'float32'})
	f['GeometryData'] = uproot.newtree({'EventId':'int32','FiberId':'int32','FiberTheta':'float32','FiberPhi':'float32'})

	f['SiPMData']['Integral'].newbasket(integral)
	f['SiPMData']['Peak'].newbasket(peak)
	f['SiPMData']['ToA'].newbasket(tstart)
	f['SiPMData']['ToP'].newbasket(ptime)
	f['SiPMData']['ToT'].newbasket(tover)

	f['GeometryData']['EventId'].newbasket(event)
	f['GeometryData']['FiberId'].newbasket(fiber)
	f['GeometryData']['FiberTheta'].newbasket(theta)
	f['GeometryData']['FiberPhi'].newbasket(phi)

def SaveWaves(fname,output):
	signals = np.vstack(output[:,4])

	sipmsettings = [size,
					cellsize,
					dcr,
					xt,
					ap,
					cellrecovery*sampling,
					trise*sampling,
					tfall*sampling,
					-20*np.log10(SNR**2),
					ccgv]

	with h5py.File('waveforms.hdf5','a') as hf:
		fname = fname+'0'
		groups = list(hf.keys())
		if fname in groups:
			i = int(groups[-1][-1])
			fname = fname[:-1]+str(i+1)
		grp = hf.create_group(fname)
		dset1 = grp.create_dataset('Waveforms',shape=(signals.shape),dtype='f',compression='lzf',chunks=(1,signals.shape[1]))
		dset2 = grp.create_dataset('SiPMSettings',shape=(len(sipmsettings),),dtype='f',compression='gzip',compression_opts=9)
		dset1[...] = signals
		dset2[...] = sipmsettings
