#!/usr/bin/env python2
import MDAnalysis
import numpy as np
from numpy.linalg import norm

u = MDAnalysis.Universe('../../ionized.psf', '../output-dcd/NPT-250-pf10ps.dcd')
#print(u)
#print(list(u.atoms[:100].residues))

sel_prot = u.select_atoms("protein and backbone")
sel_memb = u.select_atoms("segid L.* and prop z > 0.0")

dist = []
f = open("separation.dat","w")
for ts in u.trajectory:
    print("Frame: %5d, Time: %8.3f ps" % (ts.frame, u.trajectory.time))
    A = sel_prot.center_of_mass()
    B = sel_memb.center_of_mass()
    dummy = A - B
    dist.append((u.trajectory.time, norm(dummy)))
    f.write(str(ts.frame) + " " + str(norm(dummy)) + "\n")
f.close()
print(dist)
