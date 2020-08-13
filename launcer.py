from main import *

# Openig file
fname = '../Data/out.txt'
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
    # if np.all(t < SIGLEN):
    TIMES.append(t)
    OTHER.append((*L[:6],))
    if int(L[0]) % 100 == 0 and int(L[0]) > 0 and temp != L[0]:
        temp = L[0]
        print(f'Reading event: {int(L[0]):d} / {int(lines[-1].split()[0]):d}',end='\r')
    # if L[0] == '100':
    #     break
del lines

NFIB = len(TIMES)
NEVT = int(L[0])
BATCHSIZE = 100000

# Setting up results list
pool = Pool(processes=nJobs, initializer=initializeRandomPool)
res = [None] * BATCHSIZE
other = [None] * NFIB
output = np.empty((NFIB, 5), dtype='float32')
if args.wavedump:
    signals = np.empty((NFIB, SIGPTS), dtype='float32')
print('\n===> Starting simulation <===\n')

# Launching simulation
Ts = time.time()
j, k = 0, 0
for i in range(NFIB):
    if j == BATCHSIZE:
        print(f'Clearing results from cache...')
        res[-1].wait()
        # Retirieving someresults to free RAM
        for r in res:
            temp = r.get()
            output[k, :] = temp[:5]
            other[k] = temp[5]
            if args.wavedump:
                signals[k, :] = temp[6]
            k += 1
            j = 0
        print(f'Signals processed:\t{i:d}')
        print(f'Events processed:\t{other[k-1][0]} / {NEVT}\n')
    res[j] = pool.apply_async(SiPM, args=(TIMES[i], OTHER[i]))
    j += 1
pool.close()
pool.join()
# Retrieving remaining results
print(f'Clearing results from cache...')
for i in range(0, j):
    temp = res[i].get()
    output[k, :] = temp[:5]
    other[k] = temp[5]
    if args.wavedump:
        signals[k, :] = temp[6]
    k += 1
print(f'Signals processed {i:d}')
print(f'Events processed: {other[k-1][0]} / {NEVT}\n')
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
