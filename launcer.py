from main import *

pool = Pool(processes = nJobs, initializer = initializeRandomPool)
res = [None] * 50

Ts = time.time()
for i in range(50):
    times = np.ones(np.random.poisson(25)) + 30
    other = (i)
    # SiPM(np.array(times),other)
    res[i] = pool.apply_async(SiPM, args=(times, other))
pool.close()
pool.join()
Te = time.time()
print(f'Execution time: {50/(Te-Ts):.2f}s')


output = np.empty((len(res), 5))
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
