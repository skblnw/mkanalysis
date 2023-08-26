import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import matplotlib as mpl

# #############################################################################
# Generate the data

PDB='/data1/kevin/Ups1-Mdm35/Gromacs/Ups1_gmx_t1/pdb2gmx/prot.pdb'
XTC='/data1/kevin/Ups1-Mdm35/Gromacs/Ups1_gmx_t1/xtc/prot-1-8-pf10ps-1006ns-fit-50000.xtc'
traj = md.load(XTC, top=PDB)
traj.superpose(traj, 0)

top = traj.topology
sel = top.select('name CA')
sel_traj = traj.xyz[:, sel, :]

X = sel_traj.reshape(traj.n_frames, sel.shape[0] * 3)
X -= X.mean(axis=0)
X /= X.std(axis=0)

emp_cov = np.dot(X.T, X)/traj.n_frames

emp_cov_residue=np.zeros((sel.shape[0],sel.shape[0]))
for ii in range(0, sel.shape[0]):
    for jj in range(0, sel.shape[0]):
        emp_cov_residue[ii][jj] = ( emp_cov[ii*3][jj*3] + emp_cov[ii*3+1][jj*3+1] + emp_cov[ii*3+2][jj*3+2] ) / 3


# #############################################################################
# Plot the results
fig, ax = plt.subplots()
im = ax.imshow(emp_cov_residue, interpolation='gaussian', vmin=-1, vmax=1,cmap=plt.cm.RdBu_r)
ax.set_xlim(1, 169)
ax.set_ylim(1, 169)
ax.set_xticks(range(1,169,50))
ax.set_yticks(range(1,169,50))

plt.colorbar(im)


plt.tight_layout()
plt.savefig("covar.pdf")
