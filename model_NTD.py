from modeller import * 
from modeller.automodel import * 
log.verbose()    # request verbose output
env = environ()  # create a new MODELLER environment to build this model in
env.io.atom_files_directory = './:../atom_files'
a = automodel(env,
              alnfile  = 'alignm_NTD.pir',
              knowns   = ('c4onsB_','c2g57A_','PoingAA', ),
              sequence = 'query')
a.starting_model= 1
a.ending_model  = 50
a.make()
