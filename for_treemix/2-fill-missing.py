# J. Melo-Ferreira 10/04/2018
##
# Fills in sites with missing information present in other populations. Outputs total table of variants, with the same sites in all pops.
##
# Usage:
# python 2-fill-missing.py angsd-file n-individuals complete-keys > output
# Run per population
##



import sys

dict_pop = {}

with open(sys.argv[1], 'r') as f:
	header = f.readline()
	for line in f:
		if line is not '':
			columns = line.split('\t')
			key = '|'.join([columns[0], columns[1]])
			freq = float(columns[6])
			der = int(round(freq*int(sys.argv[2])*2))
			anc = (int(sys.argv[2])*2)-der
			dict_pop.update({key: str(der)+','+str(anc)})

f.close()

with open(sys.argv[3], 'r') as g:
	for coord in g:
		if coord is not '':
			coord = coord.rstrip('\n')
			if coord is not '':
				if coord in dict_pop:
					print coord+'\t'+dict_pop[coord]			 
				else:
					print coord+'\tx'

g.close()

sys.exit()

