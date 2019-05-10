import parmed as pmd
gmx = pmd.load_file('topol.top', xyz='ionized.pdb')
gmx.save('ionized.psf')
