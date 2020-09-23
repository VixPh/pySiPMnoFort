from main import *

# Openig file
fname = '../Data/BandX0/br2_seed1_timefinal.digi'
f = open(fname)
print(f'Opening file: {fname}')
lines = f.readlines()
f.close()

# Reading txt file
TIMES = []
OTHER = []
temp = 0
for line in lines:
    if not line.strip():
        continue
    L = line.split()
    t = np.array(L[6:], dtype='float32')
    TIMES.append(t)
    OTHER.append((np.float32(L[0]), np.float32(L[1] == 'Scin'), np.float32(L[2]), np.float32(L[3]), np.float32(L[4]), np.float32(L[5])))
    if int(L[0]) % 100 == 0 and int(L[0]) > 0 and temp != L[0]:
        temp = L[0]
        print(f'Reading event: {int(L[0]):d} / {int(lines[-1].split()[0]):d}', end='\r')
    if L[0] == '300':
        break
del lines

OTHER = np.array(OTHER)
NFIB = len(TIMES)
NEVT = int(L[0])
BATCHSIZE = 250000

# Setting up results list
pool = Pool(processes=nJobs, initializer=initializeRandomPool)
res = [None] * BATCHSIZE
output = np.empty(shape=(NFIB, 5), dtype='float32')
other = np.empty(shape=(NFIB, 6), dtype='float32')
if args.wavedump:
    signals = np.memmap("tmp", shape=(NFIB, SIGPTS), dtype='float32', mode = "w+")
print('\n===> Starting simulation <===\n')

# Launching simulation
Ts = time.time()
j, k = 0, 0
for i in range(NFIB):
    if j == BATCHSIZE:
        print(f'Clearing results from cache...')
        res[-1].wait()
        # Retirieving some results to free RAM
        for r in res:
            temp = r.get()
            output[k, :] = temp[:5]
            other[k, :] = temp[5]
            if args.wavedump:
                signals[k, :] = temp[6]
            k += 1
            j = 0
        print(f'Signals processed:\t{i:d}')
        print(f'Events processed:\t{int(other[k-1,0]):d} / {NEVT}\n')
    res[j] = pool.apply_async(SiPM, args=(TIMES[i], OTHER[i, :]))
    j += 1
pool.close()
pool.join()
# Retrieving remaining results
print(f'Clearing results from cache...')
for i in range(0, j):
    temp = res[i].get()
    output[k, :] = temp[:5]
    other[k,:] = temp[5]
    if args.wavedump:
        signals[k, :] = temp[6]
    k += 1
print(f'Signals processed:\t{i:d}')
print(f'Events processed:\t{int(other[k-1,0]):d} / {NEVT}\n')
Te = time.time()

print('\n===> Simulation finished <===\n')
print(f'Execution time: {(Te-Ts):.2f}s')
print(f'Events per second: {NEVT/(Te-Ts):.2f}')

integral = output[:, 0]
peak = output[:, 1]
tstart = output[:, 2]
tover = output[:, 3]
ptime = output[:, 4]

if args.graphics:
    print('\n===> Generating plots <===\n')
    somestats(output)

if args.write:
    print('\n===> Writing results <===\n')
    SaveFile(args.write, output, other)

if args.wavedump:
    print('\n===> Writing waveforms <===\n')
    SaveWaves(args.wavedump, signals)
    if(os.path.exists("temp")):
        os.remove("temp")
