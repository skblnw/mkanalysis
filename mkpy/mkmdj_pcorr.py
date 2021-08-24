#!/home/kevin/anaconda3/envs/mdtraj/bin/python
import argparse
from argparse import RawDescriptionHelpFormatter
import mdtraj as mdj
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
    'xtick.major.size': 8,
    'xtick.major.width': 3,
    'xtick.direction': 'in',
    'ytick.labelsize': 16,
    'ytick.major.size': 8,
    'ytick.major.width': 3,
    'ytick.direction': 'in',
    'errorbar.capsize': 20,
    'figure.figsize': [4, 4],
    'figure.dpi': 300
   }
mpl.rcParams.update(params)

# /----------------------------------/
# /         Argument Parser          /
# /----------------------------------/

ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
ap.add_argument('pdb', type=str, help='Structure File')
ap.add_argument('trj', type=str, help='Trajectory File')
ap.add_argument('ref', type=str, help='Reference Structure File', nargs='?', default='noref')
ap.add_argument('--res', type=int, help='Residue Starting ID', nargs='?', default=1)
ap.add_argument('--sel', type=str, help='Selection for alignment', nargs='?', default='name CA')
ap.add_argument('--xlabel', type=str, help='X-axis label', nargs='?', default=' ')
ap.add_argument('--ylabel', type=str, help='Y-axis label', nargs='?', default=' ')
ap.add_argument('--title', type=str, help='Title label', nargs='?', default=' ')
io_group = ap.add_mutually_exclusive_group(required=True)
io_group.add_argument('-o', '--output', type=str, help='PDF output file')
io_group.add_argument('-i', '--interactive', action='store_true', 
                    help='Launches an interactive matplotlib session')
parser = ap.parse_args()

# /----------------------------------/
# /            Read Data             /
# /----------------------------------/

resid_start = parser.res

if parser.ref == 'noref':
    ref = mdj.load(parser.pdb)
else:
    ref = mdj.load(parser.ref)

print("mdj> Loading trajectory")
traj = mdj.load(parser.trj,
             top=parser.pdb)
print("mdj> Finished Loading")

print("mdj> Aligning trajectory according to selection")
sel = traj.topology.select(parser.sel)
traj.superpose(ref, 0, atom_indices=sel)
print("mdj> Finished Alignment")

sel_traj = traj.xyz[:, sel, :]
sel_coor = sel_traj.reshape(traj.n_frames, sel.shape[0] * 3)

# /----------------------------------------------------/
# /               Calculations start here              /
# /----------------------------------------------------/

X = np.true_divide(np.subtract(sel_coor, np.average(sel_coor, axis=0)), sel_coor.std(axis=0))
# X = np.subtract(sel_coor, np.average(sel_coor, axis=0))
emp_cov = np.dot(X.T, X) / traj.n_frames

emp_cov_residue=np.zeros((sel.shape[0],sel.shape[0]))
for ii in range(0, sel.shape[0]):
    for jj in range(0, sel.shape[0]):
        emp_cov_residue[ii][jj] = ( emp_cov[ii*3][jj*3] + emp_cov[ii*3+1][jj*3+1] + emp_cov[ii*3+2][jj*3+2] ) / 3
        # emp_cov_residue[ii][jj] = abs(emp_cov_residue[ii][jj])


# /----------------------------------------------------/
# /                     Plotting Area                  /
# /----------------------------------------------------/

fig, ax = plt.subplots()
im = ax.imshow(emp_cov_residue, interpolation='gaussian', cmap=plt.cm.RdBu, 
                                vmin=-1, vmax=1, norm=colors.CenteredNorm(halfrange=.6),
                                origin='lower', 
                                extent=[resid_start,resid_start-1+len(sel),resid_start,resid_start-1+len(sel)]
                                )
# ax.set_xlim(430, 515)
# ax.set_ylim(430, 515)
# ax.set_xticks(range(445,515,15))
# ax.set_yticks(range(445,515,15))

cbar = plt.colorbar(im,
                    ticks=[-1,-.5,0,.5,1],
                    fraction=0.046, pad=0.04)
cbar.outline.set_linewidth(1)
cbar.ax.tick_params(length=2, width=1)

plt.tight_layout()
if parser.interactive:
    plt.show()
else:
    plt.savefig(parser.output)
