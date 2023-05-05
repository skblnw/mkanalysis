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

def calc_rmsd_align_ruvc(ulist,ref,outdir):
    try:
        os.mkdir(outdir)
    except OSError:
        print ("mda>Directory %s already exists. File will be rewritten. Be careful!" % outdir)
    else:
        print ("mda>Successfully created the directory %s " % outdir)

    ALL='name CA'
    REC1='name CA and resid 24-338'
    REC2='name CA and resid 339-590'
    WED='name CA and (resid 1-24 or resid 591-661 or resid 762-891)'
    PI='name CA and resid 662-761'
    RUVC='name CA and (resid 892-1086 or resid 1250-1300)'
    NUC='name CA and resid 1087-1249'

    rmsdlist = []
    for ii, u in enumerate(ulist):
        R = rms.RMSD(u,
                    ref,
                    select='name CA and (resid 892-1086 or resid 1250-1300)',
                    groupselections=[ALL, REC1, REC2, WED, PI, RUVC, NUC],
                    ref_frame=0).run()
        rmsdlist.append(R.rmsd[:,2])
        np.savetxt("{}/{}.dat".format(outdir,ii),np.c_[R.rmsd[:,0],R.rmsd[:,2:10]])

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
    #calc_rmsd(ulist,ref_close,outdir)
    outdir='rmsd_close_align_ruvc'
    #calc_rmsd_align_ruvc(ulist,ref_close,outdir)

    file_list = [
        '4.xtc',
        '5.xtc',
        '6.xtc',
        '8-1.xtc',
        '8-2.xtc',
        '6i1k-helixlid-DA-ca.xtc']

    ulist = [ mda.Universe('open_ca.gro',trajectory) for trajectory in file_list]
    outdir='rmsd_open'
    #calc_rmsd(ulist,ref_open,outdir)
    outdir='rmsd_open_align_ruvc'
    calc_rmsd_align_ruvc(ulist,ref_open,outdir)
