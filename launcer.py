from main import *
from tqdm import tqdm

fname = '/tmp/out.txt'
f = open(fname)
lines = f.readlines()

TIMES = []
OTHER = []
print('Reading input file')
pbar = tqdm(total=len(lines))
for l in lines:
    if not l.strip():
        continue
    L = l.split()
    TIMES.append(np.array(L[6:], dtype='float32'))
    OTHER.append((*L[:6],))
    pbar.update(1)
pbar.close()
NEVTS = int(L[0])

pool = Pool(processes=nJobs, initializer=initializeRandomPool)
res = [None] * NEVTS

print('Starting simulation')
Ts = time.time()
for i in range(NEVTS):
    res[i] = pool.apply_async(SiPM, args=(TIMES[i], OTHER[i]))
pool.close()
pool.join()
Te = time.time()
print('Simulation finished')
print(f'Execution time: {(Te-Ts):.2f}s')

output = np.empty((len(res), 5), dtype='float32')
if args.wavedump:
    signals = np.empty((len(res), SIGPTS), dtype='float32')

other = []
for i, r in enumerate(res):
    temp = r.get()
    output[i, :] = temp[:5]
    other.append(temp[5])
    if args.wavedump:
        signals[i, :] = temp[6]

integral = output[:, 0]
peak = output[:, 1]
tstart = output[:, 2]
tover = output[:, 3]
ptime = output[:, 4]

if args.graphics:
    somestats(output)

if args.write:
    SaveFile(args.write, output)

if args.wavedump:
    SaveWaves(args.wavedump, signals)
