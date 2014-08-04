### 	ATTACH MODELLED REGION TO A KNOWN REGION OF A PROTEIN 	###
## 	MUST HAVE PYTHON AND MODELLER INSTALLED			 ##
# 	change these variables for your specific project v	  #
#########################################################################
BCAT_seq = 'QDDAELATRAIPELTKLLNDEDQVVVNKAAVMVHQLSKKEASRHAIMRSPQMVSAIVRTMQNTNDVETARCTAGTLHNLSHHREGLLAIFKSGGIPALVKMLGSPVDSVLFYAITTLHNLLLHQEGAKMAVRLAGGLQKMVALLNKTNVKFLAITTDCLQILAYGNQESKLIILASGGPQALVNIMRTYTYEKLLWTTSRVLKVLSVCSSNKPAIVEAGGMQALGLHLTDPSQRLVQNCLWTLRNLSDAATKQEGMEGLLGTLVQLLGSDDINVVTCAAGILSNLTCNNYKNKMMVCQVGGIEALVRTVLRAGDREDITEPAICALRHLTSRHQEAEMAQNAVRLHYGLPVVVKLLHPPSHWPLIKATVGLIRNLALCPANHAPLREQGAIPRLVQLLVRAHQDTQRRTSMGGTQQQFVEGVRMEEIVEGCTGALHILARDVHNRIVIRGLNTIPLFVQLLYSPIENIQRVAAGVLCELAQDKEAAEAIEAEGATAPLTELLHSRNEGVATYAAAVLFRMSEDKPQDYKKRLSVELTSSLFRTEPMAWNETADLGLDIGAQGEALGYRQDDPSYRSFHSGGYGQDALGMDPMMEHEMGGHHPGADYPVDGLPDLGHAQDLMDGLPPGDSNQLAWFDTDL'
### Sequence of known model through appended region, will be the sequence of the final model
ctd_loc = "ctd-models/" 
### relative location of models to be appended, will search for all '.pdb' files in this folder
known_mod = '2Z6H.pdb'
##########################################################################


import os
from modeller import *
from modeller.automodel import *
import modeller.salign
from modeller.scripts import complete_pdb

if not os.path.isdir('alignments/'):	# make new folders if necessary
   os.makedirs('alignments/')
if not os.path.isdir('models/'):
   os.makedirs('models/')

env = environ()
env.edat.radii_factor = 0.95	# if having clashing problems increase this number at risk of bad angles
env.edat.dynamic_sphere= True
env.libs.topology.read('${LIB}/top_heav.lib')
env.libs.parameters.read('${LIB}/par.lib')		#load parameters
env.io.atom_files_directory = ['.', '${LIB}/atom_files/']
log.none()

x = 0 # a counter that will be used for naming the models later

for mod in os.listdir(ctd_loc):
    if mod.endswith(".pdb"):
       	mdl_n = complete_pdb(env, ctd_loc+mod)	
	mdl_n.write(file='CTD.pdb')		# write a temporary pdb file to be used for modelling and aligment
	
	mdl_ori = model(env)
	aln_ori = alignment(env)
	mdl_ori.build_sequence(BCAT_seq)
	aln_ori.append_sequence(BCAT_seq)	# original model with sequence of beta-catenin
	aln_ori[0].code = 'BCAT'
	aln_ori.write(file='alignments/BCAT.ali') 

	aln_final = alignment(env, file='alignments/BCAT.ali')
	mdl_arm = model(env)
	mdl_arm.read(file= known_mod)
	aln_final.append_model(mdl_arm, align_codes=known_mod, atom_files=known_mod)
	aln_final.salign(
		  gap_function=False,
		  feature_weights=(8., 0., 0., 0., 0., 0.),	# aligning the known region with the whole sequence
		  similarity_flag=True) 

		# use the newly made pdb file as a template for the whole sequence
	aln2 = alignment(env, file='alignments/BCAT.ali')
	mdl2 = model(env)
	mdl2.read(file='CTD', model_segment=('FIRST:@','END:@'))
	aln2.append_model(mdl2, align_codes='CTD')
	aln2.align()
	aln2.write(file='alignments/BCAT_CTD.ali')

		# append the CTD to the whole sequence alignment file
	aln_final.append(file='alignments/BCAT_CTD.ali', align_codes='CTD')
	aln_final.write(file='alignments/BCAT_main.ali')
	aln_final.write(file='alignments/BCAT_main.pap', alignment_format='PAP', alignment_features='INDICES CONSERVATION')
		# make pdb models for the arm region and connected CTD
	mdl_complete = automodel(env, alnfile  = 'alignments/BCAT_main.ali', knowns = (known_mod,'CTD'), sequence = 'BCAT', assess_methods=(assess.DOPE, assess.GA341, assess.normalized_dope))
	mdl_complete.starting_model = 1
 	mdl_complete.ending_model = 1 # if you increase number of models make sure you adjust the renaming convention below
	mdl_complete.make()
	x = x + 1
	os.rename('BCAT.B99990001.pdb', 'models/BCAT.B9999'+str(x)+'.pdb')	#will place models in a new folder with 
	os.rename('BCAT.V99990001', 'models/BCAT.V9999'+str(x)+'.pdb')		# their appropriate violation file
	os.rename('BCAT.D00000001', 'models/BCAT.D0000'+str(x)+'.pdb')		# and MD records

