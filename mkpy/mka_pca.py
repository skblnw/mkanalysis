from __future__ import print_function
import mdtraj as md
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

traj = md.load('/home/PHARMACY/chan.773/storage/cpf1/run_complex/md-1-100-pf100ps-fit.xtc', top='/home/PHARMACY/chan.773/Dropbox/QWD/CRISPR-Cpf1/cpf1/run_rna/output/initial.pdb')
pca1 = PCA(n_components=2)
traj.superpose(traj, 0)
reduced_cartesian = pca1.fit_transform(traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3))

plt.figure()
plt.scatter(reduced_cartesian[:, 0], reduced_cartesian[:,1], marker='x', c=traj.time)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('Cartesian coordinate PCA')
cbar = plt.colorbar()
cbar.set_label('Time [ps]')

plt.show()