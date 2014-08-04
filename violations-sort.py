## simple parse of violations files produced by automodel used to sort the best models
## produces two text files with columns of models and violations 
### MUST HAVE 'VIOLATIONS' PROFILES PRODUCED BY MODELLER AUTOMODEL

import os
h = open('violations.txt','w')
for mod in os.listdir('.'):
	if mod.startswith('query.V'): #CHANGE FOR YOUR NAMING CONVENTION
		m_file = open(mod, 'r')	
		for line in m_file:
			try:
				if line.split()[2] == 'sum':  # LOOKS FOR ROW WITH THE SUM OF VIOLATIONS !ONLY  USE APROPRIATE 'VIOLATION FILES'
					viol = line.split()[len(line.split())-1]
			except (IndexError):
				monkeybait = 0
		m_file.close()
		h.write(mod+'\t'+viol+'\n')
h.close()
p = open('violations.txt','r')		
data = p.readlines()
data.sort(key=lambda l: float(l.split()[1]))
r = open('violationsort.txt','w')
r.writelines(data)
