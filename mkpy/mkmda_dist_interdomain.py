#!/usr/bin/env python2
import MDAnalysis
import numpy as np
from numpy.linalg import norm

u = MDAnalysis.Universe('../../ionized.psf', '../output-dcd/NPT-250-pf10ps.dcd')
print(u)
#print(list(u.atoms[:100].residues))

sel = u.select_atoms("protein and backbone")
print(sel.center_of_mass())
print(sel.center_of_geometry())

sel_PH1 = u.select_atoms("segid P1 and resid 251-360")
sel_PH2 = u.select_atoms("segid P2 and resid 251-360")

dist = []
f = open("separation.dat","w")
for ts in u.trajectory[0:1000]:
    print("Frame: %5d, Time: %8.3f ps" % (ts.frame, u.trajectory.time))
    #print("Radius of Gy: %g A" % (u.atoms.radius_of_gyration()))
    PH1 = sel_PH1.center_of_mass()
    PH2 = sel_PH2.center_of_mass()
    dummy = PH1 - PH2
    dist.append((u.trajectory.time, norm(dummy)))
    f.write(str(ts.frame) + " " + str(norm(dummy)) + "\n")
f.close()
print(dist)
