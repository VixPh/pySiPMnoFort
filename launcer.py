from main import *

f = open('event.txt','r')
eventID = []
fiberID = []
theta = []
phi = []
times = []
# for l in f.readlines():
#     ls=l.strip().split(',')
#     eventID.append(int(ls[0]))
#     fiberID.append(int(ls[1]))
#     theta.append(float(ls[2]))
#     phi.append(float(ls[3]))
#     ti = np.asarray(ls[3:]).astype('float32')
#     times.append(ti)

pool = Pool(processes = nJobs ,initializer = initializeRandomPool)
t = time.time()
nEvt= 50000
for i in range(nEvt):
    #SiPM(np.array([30]*25),())
    times = np.ones(np.random.poisson(200))*25
    # other = np.asarray((eventID[i],fiberID[i],theta[i],phi[i]))
    pool.apply_async(SiPM,args=(times,()),callback=Callback)
pool.close()
pool.join()
T = time.time()
print(f"Events per second: {nEvt/(T-t):.2f}")
output = np.asarray(output)


if args.graphics:
    somestats(output)

if args.write is not None:
    SaveFile(args.write,output)

if args.wavedump is not None:
    SaveWaves(args.wavedump,output)
