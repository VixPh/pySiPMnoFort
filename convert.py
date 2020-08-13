import numpy as np
import matplotlib.pyplot as plt

f = open('br0_seed0_timefinal.txt')
lines = f.readlines()
f.close()
fout = open('out.txt','w')


for i,l in enumerate(lines):
    if l.startswith('->Event'):
        enum = int(l.split()[-1])
        continue
    if l.startswith('->Cher') or l.startswith('->Scin'):
        pass
        fout.write('\r' + str(enum) + ' ' + ' '.join(l[2:].split()[:-1]) + ' ')
    else:
        fout.write(l.split()[0] + ' ')

plt.show()
fout.close()
