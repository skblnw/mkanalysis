import matplotlib.pyplot as plt
import numpy as np
import mdshare
import pyemma

#pdb = mdshare.fetch('alanine-dipeptide-nowater.pdb', working_directory='data')
#files = mdshare.fetch('alanine-dipeptide-*-250ns-nowater.xtc', working_directory='data')
pdb = '/data2/kevin/cpf1/run_complex/pdb2gmx/prot.pdb'
files = '/data2/kevin/cpf1/run_complex/xtc/md-1-100-pf100ps-fit.xtc'
print(pdb)
print(files)

feat = pyemma.coordinates.featurizer(pdb)
feat.add_selection(feat.select_Ca())
#feat.add_selection(feat.select_Heavy())
#feat.add_backbone_torsions(periodic=False)
#feat.add_sidechain_torsions(which='chi1', cossin=True, periodic=False)
#feat.add_distances_ca(periodic=False)
# Feature - pairs
#pairs = feat.pairs(feat.select_Heavy())
#feat.add_distances(pairs, periodic=False)

data = pyemma.coordinates.load(files, features=feat)

print('We have {} features.'.format(feat.dimension()))

pca = pyemma.coordinates.pca(data, dim=2)
tica = pyemma.coordinates.tica(data, lag=3, dim=2)

pca_concatenated = np.concatenate(pca.get_output())
tica_concatenated = np.concatenate(tica.get_output())

cls_pca = pyemma.coordinates.cluster_kmeans(pca, k=100, max_iter=50, stride=10)
cls_tica = pyemma.coordinates.cluster_kmeans(tica, k=100, max_iter=50, stride=10)

its_pca = pyemma.msm.its(
    cls_pca.dtrajs, lags=[1, 2, 5, 10, 20, 50], nits=4, errors='bayes')
its_tica = pyemma.msm.its(
    cls_tica.dtrajs, lags=[1, 2, 5, 10, 20, 50], nits=4, errors='bayes')

fig, axes = plt.subplots(2, 3, figsize=(12, 6))
pyemma.plots.plot_feature_histograms(pca_concatenated, ax=axes[0, 0])
pyemma.plots.plot_feature_histograms(tica_concatenated, ax=axes[1, 0])
axes[0, 0].set_title('PCA')
axes[1, 0].set_title('TICA')
pyemma.plots.plot_density(*pca_concatenated.T, ax=axes[0, 1], cbar=False, alpha=0.1)
axes[0, 1].scatter(*cls_pca.clustercenters.T, s=15, c='C1')
axes[0, 1].set_xlabel('PC 1')
axes[0, 1].set_ylabel('PC 2')
pyemma.plots.plot_density(*tica_concatenated.T, ax=axes[1, 1], cbar=False, alpha=0.1)
axes[1, 1].scatter(*cls_tica.clustercenters.T, s=15, c='C1')
axes[1, 1].set_xlabel('IC 1')
axes[1, 1].set_ylabel('IC 2')
pyemma.plots.plot_implied_timescales(its_pca, ax=axes[0, 2], units='ps')
pyemma.plots.plot_implied_timescales(its_tica, ax=axes[1, 2], units='ps')
axes[0, 2].set_ylim(1, 2000)
axes[1, 2].set_ylim(1, 2000)
fig.tight_layout()

plt.savefig('msm.png')
