from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()

env = environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

a = MyModel(env,
				alnfile  = 'alignment.seq',
				knowns   = 'chainA',
				sequence = 'chainA_fill')
a.starting_model= 1
a.ending_model  = 1

a.make()
