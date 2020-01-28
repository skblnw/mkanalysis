import matplotlib.pyplot as plt
import numpy as np
import mdshare
import pyemma

#pdb = mdshare.fetch('alanine-dipeptide-nowater.pdb', working_directory='data')
#files = mdshare.fetch('alanine-dipeptide-*-250ns-nowater.xtc', working_directory='data')
#feat.add_selection(feat.select_Heavy())
#feat.add_sidechain_torsions(which='chi1', cossin=True, periodic=False)
#feat.add_selection(feat.select_Ca())
#feat.add_distances_ca(periodic=False)

pdb = '/data1/kevin/Ups1-Mdm35/Gromacs/Ups1_gmx_t1/pdb2gmx/prot.pdb'
files = ['/data1/kevin/Ups1-Mdm35/Gromacs/Ups1_gmx_t1/xtc/prot-1-8-pf10ps-1006ns-fit-50000.xtc', '/data1/kevin/Ups1-Mdm35/Gromacs/Ups1_gmx_t2/xtc/prot-1-8-pf10ps-915ns-fit-50000.xtc']
feat = pyemma.coordinates.featurizer(pdb)
feat.add_backbone_torsions(periodic=False)
data = pyemma.coordinates.load(files, features=feat)


print('We have {} features.'.format(feat.dimension()))

pca = pyemma.coordinates.pca(data)
pca_concatenated = np.concatenate(pca.get_output())

fig, axes = plt.subplots(1, 2, figsize=(12,5))
pyemma.plots.plot_density(*pca_concatenated[:, :2].T, ax=axes[0], cbar=False, logscale=True)
pyemma.plots.plot_free_energy(*pca_concatenated[:, :2].T, ax=axes[1], legacy=False)
for ax in axes.flat[1:]:
    ax.set_xlabel('PC 1')
    ax.set_ylabel('PC 2')
fig.tight_layout()
#plt.savefig('gmx_pca.png')

tica = pyemma.coordinates.tica(data)
tica_concatenated = np.concatenate(tica.get_output())

fig, axes = plt.subplots(1, 2, figsize=(12,5))
pyemma.plots.plot_density(*tica_concatenated[:, :2].T, ax=axes[0], cbar=False, logscale=True)
pyemma.plots.plot_free_energy(*tica_concatenated[:, :2].T, ax=axes[1], legacy=False)
for ax in axes.flat[1:]:
    ax.set_xlabel('IC 1')
    ax.set_ylabel('IC 2')
fig.tight_layout()
#plt.savefig('gmx_tica.png')

fig, ax = plt.subplots(figsize=(10,3))
ax.plot(tica_concatenated[:,1], label='TICA')
#ax.plot(pca_concatenated[:,:1], label='PCA')
ax.legend()
fig.tight_layout()
#plt.savefig('gmx_tica_ic1.png')

fig, ax = plt.subplots(figsize=(10,3))
ax.plot(tica_concatenated[:,2], label='TICA')
#ax.plot(pca_concatenated[:,2], label='PCA')
ax.legend()
fig.tight_layout()
#plt.savefig('gmx_tica_ic2.png')

plt.show()
