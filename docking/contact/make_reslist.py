import os
import numpy as np
A=np.zeros(1307)
f=open("namelist", "r")
for ii in f:
    dir=os.path.dirname(ii)
    filename=os.path.join(dir, "out.contact")

    f=open(filename,"r")
    for ii in f:
        ii=ii.split()
        jj=list(map(int, ii))
    for kk in jj:
        A[kk]=A[kk]+1

print(A)
np.savetxt("reslist.txt", A)
