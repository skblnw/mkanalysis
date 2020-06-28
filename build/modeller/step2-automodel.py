from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()

env = environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = automodel(env,
    alnfile  = 'alignment.seq',
    knowns   = 'prot',
    sequence = 'prot_fill',
    assess_methods=(assess.DOPE,
                    assess.GA341))
a.starting_model= 1
a.ending_model  = 1

a.make()
