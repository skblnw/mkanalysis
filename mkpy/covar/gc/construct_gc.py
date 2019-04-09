import numpy as np
import matplotlib.pyplot as plt

f = open('correl.dat', 'r')
firstline = f.readline()
size=int(firstline.split()[0])

data_raw = [line.strip().split() for line in f.readlines()]
data_raw[-1] = data_raw[-1][:-1]
data = [data for xx in data_raw for data in xx]

npa = np.array(data, dtype='float')
new = npa.reshape(size,size)
#np.savetxt(new, 'out')

fig, ax = plt.subplots()
plt.imshow(new, interpolation='gaussian', vmin=0, vmax=1,
           cmap=plt.cm.RdBu_r)
ax.set_xlim(1,169)
ax.set_ylim(1,169)
plt.colorbar()
plt.show()
