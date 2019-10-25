# J. Melo-Ferreira 10/04/2018
##
# Opens angsd file and outputs chr|site key
##
# Usage:
# python 0-keys.py angsd-file > output 
# Run per population


import sys

keys = []


with open(sys.argv[1], 'r') as f:
	header = f.readline()
	for line in f:
		if line is not '':
			columns = line.split('\t')
			key = '|'.join([columns[0], columns[1]])
			print key
		else:
			pass

f.close()

sys.exit()



