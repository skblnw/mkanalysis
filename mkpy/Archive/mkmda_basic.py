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
        np.savetxt("{}/{}.dat".format(outdir,ii),np.c_[R.rmsd[:,0],R.rmsd[:,2]])


# MAIN

if __name__ == '__main__':

    ref_close = mda.Universe('close_ca.gro')
    ref_open = mda.Universe('open_ca.gro')


    file_list = ['1.xtc',
        '2.xtc',
        '6gtd-closelid-BD-4-ca.xtc',
        '6gtd-closelid-BD-5-ca.xtc',
        '6gtd-closelid-BD-6-ca.xtc',
        '6gtd-closelid-BD-7-ca.xtc',
        '6gtd-closelid-BD-9-ca.xtc',
        '6gtd-closelid-BD-10-ca.xtc',
        '6gtd-helixlid-BD-4-ca.xtc',
        '6gtd-helixlid-BD-5-ca.xtc',
        '6gtd-helixlid-BD-9-ca.xtc',
        '6gtd-helixlid-BD-11-ca.xtc',
        '3.xtc'
        ]
    ulist = [ mda.Universe('close_ca.gro',trajectory) for trajectory in file_list]
    outdir='rmsd_close'
    calc_rmsd(ulist,ref_close,outdir)

    file_list = [
        '4.xtc',
        '5.xtc',
        '6.xtc',
        '8-1.xtc',
        '8-2.xtc',
        '6i1k-helixlid-DA-ca.xtc']

    ulist = [ mda.Universe('open_ca.gro',trajectory) for trajectory in file_list]
    outdir='rmsd_open'
    calc_rmsd(ulist,ref_open,outdir)
