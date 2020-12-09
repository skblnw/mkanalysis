import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align,rms,pca,contacts
import os


def calc_native_contact(ulist,ref,outdir):
    try:
        os.mkdir(outdir)
    except OSError:
        print ("mda> Directory %s already exists. File will be rewritten. Be careful!" % outdir)
    else:
        print ("mda> Successfully created the directory %s " % outdir)



    for ii, u in enumerate(ulist):
        sel = "name CA and resid 1008-1022"
        refg = ref.select_atoms(sel)

        print("mda> Aligning trajectory according to selection")
        aligner = align.AlignTraj(u,
                                  ref,
                                  select='name CA and resid 1008-1022',
                                  in_memory=True).run()

        print("Computing native contacts")
        ca = contacts.Contacts(u,
                                selection=(sel, sel),
                                refgroup=(refg, refg),
                                radius=4.5,
                                method='soft_cut',
                                kwargs={'beta': 5.0,'lambda_constant': 1.5}
                                ).run()
        np.savetxt("{0}/{1}.dat".format(outdir,ii),ca.timeseries)

def calc_rmsd_lid(ulist,ref,outdir):
    try:
        os.mkdir(outdir)
    except OSError:
        print ("mda>Directory %s already exists. File will be rewritten. Be careful!" % outdir)
    else:
        print ("mda>Successfully created the directory %s " % outdir)

    rmsdlist = []
    for ii, u in enumerate(ulist):
        R = rms.RMSD(u,
                    ref,
                    select='name CA and resid 1008-1022',
                    ref_frame=0).run()
        rmsdlist.append(R.rmsd[:,2])
        np.savetxt("{}/{}.dat".format(outdir,ii),np.c_[R.rmsd[:,0],R.rmsd[:,2]])

# MAIN

if __name__ == '__main__':

    ref_close = mda.Universe('close_ca.gro')
    ref_open = mda.Universe('open_ca.gro')
    ref_helixlid = mda.Universe('helixlid.pdb')

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
                '3.xtc',
                '4.xtc',
                '5.xtc',
                '6.xtc',
                '8-1.xtc',
                '8-2.xtc',
                '6i1k-helixlid-DA-ca.xtc']
    ulist = [ mda.Universe('close_ca.gro',trajectory) for trajectory in file_list]
    outdir='nativecontact_closelid'
    # calc_rmsd_lid(ulist,ref_helixlid,outdir)

    file_list = ['pdb/6gtc.gro',
                'pdb/6gtd.gro',
                'pdb/6gte.gro',
                'pdb/6i1k_closelid.gro',
                'pdb/6i1k_helixlid.gro']
    ulist = [ mda.Universe(topology) for topology in file_list]
    outdir='nativecontact_closelid_pdb'
    calc_rmsd_lid(ulist,ref_helixlid,outdir)
