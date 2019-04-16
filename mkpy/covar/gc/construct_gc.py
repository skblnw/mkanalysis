import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import matplotlib as mpl

f = open('correl_apo.dat', 'r')
firstline = f.readline()
size=int(firstline.split()[0])

data_raw = [line.strip().split() for line in f.readlines()]
data_raw[-1] = data_raw[-1][:-1]
data = [data for xx in data_raw for data in xx]

npa = np.array(data, dtype='float')
new = npa.reshape(size,size)
#np.savetxt(new, 'out')

fig, ax = plt.subplots()
im = ax.imshow(new, interpolation='gaussian', vmin=0, vmax=1,cmap=plt.cm.jet)
ax.set_xlim(0, 1306)
ax.set_ylim(0, 1306)
plt.colorbar(im)

plt.tight_layout()
plt.savefig("gc_apo.pdf")
#plt.show()
