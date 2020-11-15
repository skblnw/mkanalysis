from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()

class MyModel(automodel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['K', 'L', 'M', 'N'])
env = environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = automodel(env,
    alnfile  = 'prot.seq',
    knowns   = 'prot',
    sequence = 'prot_fill',
    assess_methods=(assess.DOPE,
                    assess.GA341))
a.starting_model= 1
a.ending_model  = 3

a.make()
