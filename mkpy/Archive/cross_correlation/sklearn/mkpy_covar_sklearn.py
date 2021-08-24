from __future__ import print_function
print(__doc__)
# author: Gael Varoquaux <gael.varoquaux@inria.fr>
# License: BSD 3 clause
# Copyright: INRIA

import numpy as np
from scipy import linalg
from sklearn.datasets import make_sparse_spd_matrix
from sklearn.covariance import GraphLassoCV, ledoit_wolf
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import mdtraj as md

# #############################################################################
# Generate the data
n_samples = 5
n_features = 4

traj = md.load('/home/PHARMACY/chan.773/storage/cpf1/run_rna/md-1-80-pf100ps-fit.xtc', 
           top='/home/PHARMACY/chan.773/Dropbox/QWD/CRISPR-Cpf1/cpf1/run_rna/output/initial.pdb')
traj.superpose(traj, 0)

top = traj.topology
sel = top.select('name CA')
sel_traj = traj.xyz[:, sel, :]

X = sel_traj.reshape(traj.n_frames, sel.shape[0] * 3)
X -= X.mean(axis=0)
X /= X.std(axis=0)

pca=PCA().fit(X)
cov = pca.get_covariance()

# #############################################################################
# Plot the results
plt.figure(figsize=(10, 6))
plt.subplots_adjust(left=0.02, right=0.98)

# plot the covariances
covs = [('True', cov)]
vmax = cov.max()
for i, (name, this_cov) in enumerate(covs):
    plt.subplot(2, 4, i + 1)
    plt.imshow(this_cov, interpolation='nearest', vmin=-vmax, vmax=vmax,
               cmap=plt.cm.RdBu_r)
    plt.xticks(())
    plt.yticks(())
    plt.title('%s covariance' % name)

plt.show()
