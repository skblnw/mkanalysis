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
#feat.add_selection(feat.select_Heavy())
#feat.add_backbone_torsions(periodic=False)
#feat.add_sidechain_torsions(which='chi1', cossin=True, periodic=False)
#feat.add_selection(feat.select_Ca())
#feat.add_distances_ca(periodic=False)
#data = pyemma.coordinates.load(files, features=feat)

feat = pyemma.coordinates.featurizer(pdb)
feat.add_backbone_torsions(periodic=False)
print('We have {} features.'.format(feat.dimension()))
data = pyemma.coordinates.load(files, features=feat)
score_phi_psi = pyemma.coordinates.vamp(data[:-1], dim=1).score(
        test_data=data[:-1],
        score_method='VAMP2')

feat = pyemma.coordinates.featurizer(pdb)
feat.add_selection(feat.select_Ca())
print('We have {} features.'.format(feat.dimension()))
data = pyemma.coordinates.load(files, features=feat)
score_heavy_atoms = pyemma.coordinates.vamp(data[:-1], dim=2).score(
        test_data=data[:-1],
        score_method='VAMP2')


print('VAMP2-score backbone torsions: {:f}'.format(score_phi_psi))
print('VAMP2-score CA: {:f}'.format(score_heavy_atoms))


