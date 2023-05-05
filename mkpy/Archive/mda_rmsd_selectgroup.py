import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align,rms
import os


def calc_rmsd(ulist,ref,outdir):
    try:
        os.mkdir(outdir)
    except OSError:
        print ("mda>Directory %s already exists. File will be rewritten. Be careful!" % outdir)
    else:
        print ("mda>Successfully created the directory %s " % outdir)


    rmsdlist = []
    for ii, u in enumerate(ulist):
        R = rms.RMSD(u,ref,ref_frame=0).run()
        rmsdlist.append(R.rmsd[:,2])
        np.savetxt("{}/{}.dat".format(outdir,ii+1),np.c_[R.rmsd[:,0],R.rmsd[:,2]])

def calc_rmsd_selectgroup(ulist,ref,outdir):
    try:
        os.mkdir(outdir)
    except OSError:
        print ("mda>Directory %s already exists. File will be rewritten. Be careful!" % outdir)
    else:
        print ("mda>Successfully created the directory %s " % outdir)

    ALL='name CA'
    PROA='name CA and resid 4-197 and segid A'
    PROB='name CA and resid 4-197 and segid B'

    rmsdlist = []
    for ii, u in enumerate(ulist):
        R = rms.RMSD(u,
                    ref,
                    select='name CA and resid 4-197 and segid A',
                    groupselections=[ALL, PROA, PROB],
                    ref_frame=0).run()
        rmsdlist.append(R.rmsd[:,2])
        np.savetxt("{}/{}.dat".format(outdir,ii+1),np.c_[R.rmsd[:,0],R.rmsd[:,3:]])

# MAIN

if __name__ == '__main__':

    ref = mda.Universe('stripped.pdb')


    file_list = ['1.xtc',
        '2.xtc',
        '3.xtc',
        '4.xtc',
        '5.xtc',
        '6.xtc'
        ]
    ulist = [ mda.Universe('stripped.pdb',trajectory) for trajectory in file_list]
    outdir = 'rmsd'
    #calc_rmsd(ulist,ref_close,outdir)
    outdir = 'rmsd_separatechain'
    calc_rmsd_selectgroup(ulist,ref,outdir)

