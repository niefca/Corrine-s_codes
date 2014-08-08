from modeller import * 
from modeller.automodel import * 
log.verbose()    # request verbose output
env = environ()  # create a new MODELLER environment to build this model in
env.io.atom_files_directory = './:../atom_files'
a = automodel(env,
              alnfile  = 'newmodeller.pir', # must have this alignment file 
              knowns   = ('d1c7ka_','PoingAA', ),
              sequence = 'query') #name of models
a.starting_model= 1
a.ending_model  = 50 # number of models made
a.make()
