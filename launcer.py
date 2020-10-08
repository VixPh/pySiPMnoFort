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
    if L[0] == '10':
        break
del lines

OTHER = np.array(OTHER)
INPUT = zip(TIMES, OTHER)
NFIB = len(TIMES)
NEVT = int(L[0])
del TIMES
del OTHER

# Setting up results arrays
pool = Pool(processes=nJobs, initializer=initializeRandomPool, maxtasksperchild=16384)
output = np.empty(shape=(NFIB, 5), dtype='float32')
other = np.empty(shape=(NFIB, 6), dtype='float32')
if args.wavedump:
    # Using memmap since this array will be VERY big
    signals = np.memmap("tmp", shape=(NFIB, SIGPTS), dtype='float32', mode="w+")
print('\n===> Starting simulation <===\n')

# Launching simulation
Ts = time.time()
res = pool.starmap_async(SiPM, INPUT)
pool.close()
pool.join()
Te = time.time()
res = res.get()

for i, r in enumerate(res):
    output[i, :] = r[:5]
    other[i, :] = r[5]

print('\n===> Simulation finished <===\n')
print(f'Execution time: {(Te-Ts):.2f}s')
print(f'Events per second: {NEVT/(Te-Ts):.2f}')

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
