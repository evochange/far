# J. Melo-Ferreira 10/04/2018
##
# Filters SNPs for distance of N bp (suggestion 20000), outputs resulting table of sites and treemix input.
# Combined-files from different chromosomes can be concatenated prior to using this script.
##
# Usage:
# python 4-filter-and-output-treemix.py combined-file n-distance-to-filter out-sites-table out-treemix-input
##




import sys

chrom = []
block = int(sys.argv[2])

out_sites = open(sys.argv[3], 'a')
out_treemix = open(sys.argv[4], 'a')

with open(sys.argv[1], 'r') as f:
	header = f.readline().rstrip('\n')
	out_sites.write(header + '\n')
	h_elements = header.split('\t')
	out_treemix.write(' '.join(h_elements[2:])+'\n')	
	for line in f:
		if line is not '':
			line = line.rstrip('\n')
			element = line.split('\t')
			if element[0] not in chrom:
				chrom = [element[0]]
				out_sites.write(line+'\n')
				out_treemix.write(' '.join(element[2:])+'\n')
			else:
				if int(element[1]) > block:
					block = block+int(sys.argv[2])
					out_sites.write(line+'\n')
					out_treemix.write(' '.join(element[2:])+'\n')

out_sites.close()
out_treemix.close()
f.close()
				
sys.exit()

