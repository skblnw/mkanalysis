import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mdtraj as mdj
from sklearn.decomposition import PCA
import os


def pca_align_ruv(t, ref, n_components):
    pc = PCA(n_components=n_components)
    t.superpose(ref,atom_indices=t.topology.select('resid 892 to 1086 or resid 1250 to 1300'))
    return pc.fit(t.xyz.reshape(t.n_frames, t.n_atoms * 3))

def pca_projection(t, pc, outdir):
    try:
        os.mkdir(outdir)
    except OSError:
        print ("mda> Directory %s already exists. File will be rewritten. Be careful!" % outdir)
    else:
        print ("mda> Successfully created the directory %s " % outdir)

    for ii, t in enumerate(tlist):
        np.savetxt("{0}/{1}.dat".format(outdir,ii),pc.transform(t.xyz.reshape(t.n_frames, t.n_atoms * 3)))

def pca_align_ruv2(u,ref):

    print("mda> Aligning trajectory")
    aligner = align.AlignTraj(u,
                              ref,
                              select='name CA and (resid 892-1086 or resid 1250-1300)',
                              in_memory=True).run()
    print("mda> Doing PCA")
    pc = pca.PCA(u,
                align=False,
                mean=None).run()

    return pc

def pca_projection2(ulist,pc,outdir):
    try:
        os.mkdir(outdir)
    except OSError:
        print ("mda> Directory %s already exists. File will be rewritten. Be careful!" % outdir)
    else:
        print ("mda> Successfully created the directory %s " % outdir)

    for ii, u in enumerate(ulist):
        sel = u.select_atoms('name CA')
        np.savetxt("{0}/{1}.dat".format(outdir,ii),pc.transform(sel, n_components=1))

# MAIN

if __name__ == '__main__':

    # mdj references
    ref_close = mdj.load('close_ca.gro')

    file_list = ['1.xtc',
                '2.xtc',
                '8-2.xtc']
    t = mdj.load(file_list, top='close_ca.gro')
    n_components=2
    pc = pca_align_ruv(t,ref_close,n_components) 

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
    tlist = [ mdj.load(trajectory, top='close_ca.gro') for trajectory in file_list]
    outdir='pca_align_ruv'
    pca_projection(tlist, pc, outdir)

