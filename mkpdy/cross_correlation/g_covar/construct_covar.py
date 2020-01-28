import numpy as np
import matplotlib.pyplot as plt

f = open('covar.dat', 'r')
size=507

data_raw = [line.strip().split() for line in f.readlines()]
data = [data for xx in data_raw for data in xx]

npa = np.array(data, dtype='float')
new = npa.reshape(size,size)
#np.savetxt(new, 'out')

vmax = new.max()
plt.imshow(new, interpolation='gaussian', vmin=-1, vmax=1,
           cmap=plt.cm.RdBu_r)
plt.colorbar()
plt.show()
