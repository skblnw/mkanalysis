import argparse
from argparse import RawDescriptionHelpFormatter
import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style
import matplotlib.font_manager
from matplotlib import colors
params = {
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',
    'font.size': 16,
    'axes.labelsize': 16,
    'axes.linewidth': 2,
    'axes.xmargin': 0,
    'axes.ymargin': 0,
    'lines.linewidth' : 2,
    'legend.fontsize': 16,
    'xtick.labelsize': 16,
    'xtick.major.size': 10,
    'xtick.major.width': 3,
    'xtick.direction': 'in',
    'ytick.labelsize': 16,
    'ytick.major.size': 10,
    'ytick.major.width': 3,
    'ytick.direction': 'in',
    'errorbar.capsize': 20,
    'figure.figsize': [8, 8],
    'figure.dpi': 300
   }
mpl.rcParams.update(params)


def align_to_ref(u,ref,seltext):
    print("mda> Aligning trajectory according to selection")
    alignment = align.AlignTraj(u,ref,
                              select=seltext,
                              in_memory=True)
    alignment.run()
    return u


if __name__ == '__main__':

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    ap.add_argument('pdb', type=str, help='Structure File')
    ap.add_argument('trj', type=str, help='Trajectory File')
    ap.add_argument('ref', type=str, help='Reference Structure File', nargs='?', default='noref')
    ap.add_argument('--resid_start', type=int, help='Residue Starting ID', nargs='?', default=1)
    ap.add_argument('--xlabel', type=str, help='X-axis label', nargs='?', default=' ')
    ap.add_argument('--ylabel', type=str, help='Y-axis label', nargs='?', default=' ')
    ap.add_argument('--title', type=str, help='Title label', nargs='?', default=' ')
    io_group = ap.add_mutually_exclusive_group(required=True)
    io_group.add_argument('-o', '--output', type=str, help='PDF output file')
    io_group.add_argument('-i', '--interactive', action='store_true', 
                        help='Launches an interactive matplotlib session')
    parser = ap.parse_args()

    resid_start = parser.resid_start

    if parser.ref == 'noref':
        ref = mda.Universe(parser.pdb)
    else:
        ref = mda.Universe(parser.ref)

    u = mda.Universe(parser.pdb,parser.trj, in_memory=True)

    seltext = 'segid B and name CA'
    sel = u.select_atoms('segid B and name CA')
    u = align_to_ref(u, ref, seltext)

    traj = np.stack([ ts.positions for ts in u.trajectory ])
    sel_traj = traj[:, sel.indices, :]
    sel_coor = sel_traj.reshape(traj.shape[0], len(sel) * 3)

    X = np.true_divide(np.subtract(sel_coor, np.average(sel_coor, axis=0)), sel_coor.std(axis=0))
    # X = np.subtract(sel_coor, np.average(sel_coor, axis=0))
    emp_cov = np.dot(X.T, X) / traj.shape[0]

    emp_cov_residue=np.zeros((len(sel),len(sel)))
    for ii in range(0, len(sel)):
        for jj in range(0, len(sel)):
            emp_cov_residue[ii][jj] = ( emp_cov[ii*3][jj*3] + emp_cov[ii*3+1][jj*3+1] + emp_cov[ii*3+2][jj*3+2] ) / 3
            emp_cov_residue[ii][jj] = abs(emp_cov_residue[ii][jj])


    ##############################################################################
    # Plot the results
    fig, ax = plt.subplots()
    im = ax.imshow(emp_cov_residue, interpolation='gaussian', 
                                    vmin=0, vmax=1, cmap=plt.cm.RdBu_r, 
                                    origin='lower', 
                                    extent=[resid_start,resid_start-1+len(sel),resid_start,resid_start-1+len(sel)])
    # ax.set_xlim(1, len(sel))
    # ax.set_ylim(1, len(sel))
    # ax.set_xticks(range(1,len(sel),50))
    # ax.set_yticks(range(1,len(sel),50))

    plt.colorbar(im)
    plt.tight_layout()

    if parser.interactive:
        plt.show()
    else:
        plt.savefig(outfile)
