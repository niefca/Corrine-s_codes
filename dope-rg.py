import os
from MDAnalysis import *
from MDAnalysis.analysis.align import *
import numpy
import numpy.linalg
from modeller import *
from modeller.scripts import complete_pdb


env = environ()
env.edat.radii_factor = 0.92
env.edat.dynamic_sphere= True
env.libs.topology.read('${LIB}/top_heav.lib')
env.libs.parameters.read('${LIB}/par.lib')	#load parameters
env.io.atom_files_directory = ['.', '${LIB}/atom_files/']
log.verbose()
table = open('analysis.txt', 'w')
for x in os.listdir('.'):
	if x.endswith('.pdb'):
		new_mod = complete_pdb(env, x)
		dope = new_mod.assess_normalized_dope()

		bcat = MDAnalysis.Universe(x)
		bcat_sel = bcat.selectAtoms('name CA or backbone')	
		bcat_RG = bcat_sel.radiusOfGyration()
		
		table.write(str(x)+'\t'+str(dope)+'\t'+str(bcat_RG)+'\n')
table.close()

