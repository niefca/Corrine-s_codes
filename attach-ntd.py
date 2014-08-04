### YOUR VARIABLES ###
FULL_seq = 'MATQADLMELDMAMEPDRKAAVSHWQQQSYLDSGIHSGATTTAPSLSGKGNPEEEDVDTSQVLYEWEQGFSQSFTQEQVADIDGQYAMTRAQRVRAAMFPETLDEGMQIPSTQFDAAHPTNVQRLAEPSQMLKHAVVNLINYQDDAELATRAIPELTKLLNDEDQVVVNKAAVMVHQLSKKEASRHAIMRSPQMVSAIVRTMQNTNDVETARCTAGTLHNLSHHREGLLAIFKSGGIPALVKMLGSPVDSVLFYAITTLHNLLLHQEGAKMAVRLAGGLQKMVALLNKTNVKFLAITTDCLQILAYGNQESKLIILASGGPQALVNIMRTYTYEKLLWTTSRVLKVLSVCSSNKPAIVEAGGMQALGLHLTDPSQRLVQNCLWTLRNLSDAATKQEGMEGLLGTLVQLLGSDDINVVTCAAGILSNLTCNNYKNKMMVCQVGGIEALVRTVLRAGDREDITEPAICALRHLTSRHQEAEMAQNAVRLHYGLPVVVKLLHPPSHWPLIKATVGLIRNLALCPANHAPLREQGAIPRLVQLLVRAHQDTQRRTSMGGTQQQFVEGVRMEEIVEGCTGALHILARDVHNRIVIRGLNTIPLFVQLLYSPIENIQRVAAGVLCELAQDKEAAEAIEAEGATAPLTELLHSRNEGVATYAAAVLFRMSEDKPQDYKKRLSVELTSSLFRTEPMAWNETADLGLDIGAQGEALGYRQDDPSYRSFHSGGYGQDALGMDPMMEHEMGGHHPGADYPVDGLPDLGHAQDLMDGLPPGDSNQLAWFDTDL'
	# entire sequence of model to be produced
n_mod = '../NTD/query.B99990023.pdb'	
	# relative location and name of model to be appended
c_mod = '../attachCTD/BCAT.B99990009.pdb'
	# relative location and name of model to be appended to 
systm = 'BCAT'
	# name of the completed system, will be used to name final models
### END OF YOUR VARIABLES ###

import os
import random
from modeller import *
from modeller.automodel import *
import modeller.salign
from modeller.scripts import complete_pdb

if not os.path.isdir('alignments/'):	# make new folders if necessary
   os.makedirs('alignments/')
if not os.path.isdir('models/'):
   os.makedirs('models/')


env = environ()
env.edat.radii_factor = 0.95 #change if encountering clashing problems at risk of bad angles
env.edat.dynamic_sphere= True
env.libs.topology.read('${LIB}/top_heav.lib')
env.libs.parameters.read('${LIB}/par.lib')	#load parameters
env.io.atom_files_directory = ['.', '${LIB}/atom_files/']
log.none()
x=0


mdl_n = complete_pdb(env, n_mod)
sel = selection()
sel.add(mdl_n)
#sel.translate([41,-4,78]) ##include a tailored version if models consistently overlap
#sel.rotate_mass_center([0,1,0],45) ## fine adjustments if you think that rotation is necessary to get best stitching #[x,y,z], degree of rotation
sel.write(file='NTD.pdb', model_format='PDB') #temporary file to be used as appended region

mdl_ori = model(env)
aln_ori = alignment(env)
mdl_ori.build_sequence(FULL_seq)
aln_ori.append_sequence(FULL_seq)	# original model with sequence of beta-catenin
aln_ori[0].code = systm
aln_ori.write(file='alignments/'+systm+'.ali') 

aln_final = alignment(env, file='alignments/'+systm+'.ali')
mdl_arm = model(env)
mdl_arm.read(file=c_mod)
aln_final.append_model(mdl_arm, align_codes='ARMC', atom_files=c_mod)
aln_final.salign(
	  gap_function=False,
	  feature_weights=(8., 0., 0., 0., 0., 0.),	# aligning the arm region with the whole sequence
	  similarity_flag=True) 

	# use the newly made pdb file as a template for the whole sequence
aln2 = alignment(env, file='alignments/'+systm+'.ali')
mdl2 = model(env)
mdl2.read(file='NTD', model_segment=('FIRST:@','END:@'))
aln2.append_model(mdl2, align_codes='NTD')
aln2.align()
aln2.write(file='alignments/'+systm+'_NTD.ali')

	# append the CTD to the whole sequence alignment file
aln_final.append(file='alignments/'+systm+'_NTD.ali', align_codes='NTD')
aln_final.write(file='alignments/'+systm+'_main.ali')
aln_final.write(file='alignments/'+systm+'_main.pap', alignment_format='PAP', alignment_features='INDICES CONSERVATION')
	# make pdb models for the arm region and connected CTD
mdl_complete = automodel(env, alnfile  = 'alignments/'+systm+'_main.ali', knowns = ('ARMC','NTD'), sequence = systm)
mdl_complete.starting_model = 1
 mdl_complete.ending_model = 1 # if you increase number of models make sure you adjust the renaming convention below
mdl_complete.make()
x = x + 1
os.rename(systm+'.B99990001.pdb', 'models/'+systm+'.B9999'+str(x)+'.pdb') #will place models in a new folder with
os.rename(systm+'.V99990001', 'models/'+systm+'.V9999'+str(x)+'.pdb') # their appropriate violation file
os.rename(systm+'.D00000001', 'models/'+systm+'.D0000'+str(x)+'.pdb') # and MD records


		
