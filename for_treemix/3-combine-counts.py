# J. Melo-Ferreira 10/04/2018
##
# Combines complete file from each population into a single table.
# Sites with missing data are removed.
# Invariable sites are removed.
# Sites with same minor allele count <N (suggestion 4) in all timidus populations and =0 in americanus are removed
##
# Usage:
# python 3-combine-counts.py complete-ALPS complete-FSCA complete-FARO complete-AMER n-minor-count > output
##


import sys
import itertools
from itertools import izip

print 'CHR\tSITE\tALPS\tFSCA\tIREL\tFARO\tAME'

ALPS_file = open(sys.argv[1], 'r')
FSCA_file = open(sys.argv[2], 'r')
FARO_file = open(sys.argv[3], 'r')
AMER_file = open(sys.argv[4], 'r')

for line_ALPS, line_FSCA, line_FARO, line_AMER in izip(ALPS_file, FSCA_file, IREL_file, FARO_file, AMER_file):
	site_ALPS = line_ALPS.split('\t')[0]
	chromosome = site_ALPS.split('|')[0]
	coordinate = site_ALPS.split('|')[1]
	var_ALPS = line_ALPS.split('\t')[1].rstrip('\n')	
	site_FSCA = line_FSCA.split('\t')[0]
	var_FSCA = line_FSCA.split('\t')[1].rstrip('\n')
	site_FARO = line_FARO.split('\t')[0]
	var_FARO = line_FARO.split('\t')[1].rstrip('\n')
	site_AMER = line_AMER.split('\t')[0]
	var_AMER = line_AMER.split('\t')[1].rstrip('\n')
	variants = [var_ALPS, var_FSCA, var_FARO, var_AMER]
	if 'x' not in variants:
		der_ALPS = var_ALPS.split(',')[0]
		anc_ALPS = var_ALPS.split(',')[1]
		der_FSCA = var_FSCA.split(',')[0]
		anc_FSCA = var_FSCA.split(',')[1]
		der_FARO = var_FARO.split(',')[0]
		anc_FARO = var_FARO.split(',')[1]
		der_AMER = var_AMER.split(',')[0]
		anc_AMER = var_AMER.split(',')[1]
		varlist_ALPS = [int(der_ALPS), int(anc_ALPS)]
		varlist_FSCA = [int(der_FSCA), int(anc_FSCA)]
		varlist_FARO = [int(der_FARO), int(anc_FARO)]
		varlist_AMER = [int(der_AMER), int(anc_AMER)]
		if site_ALPS == site_FSCA == site_FARO == site_AMER:
			if der_ALPS == der_FSCA == der_FARO == der_AMER == '0':
				pass
			else:
				if  anc_ALPS == anc_FSCA == anc_FARO == anc_AMER == '0':	
					pass
				else:
					if min(int(der_ALPS), int(anc_ALPS))>=int(sys.argv[5]):
						print chromosome+'\t'+coordinate+'\t'+var_ALPS+'\t'+var_FSCA+'\t'+var_FARO+'\t'+var_AMER
					elif min(int(der_FSCA), int(anc_FSCA))>=int(sys.argv[5]):
                                                print chromosome+'\t'+coordinate+'\t'+var_ALPS+'\t'+var_FSCA+'\t'+var_FARO+'\t'+var_AMER
					elif min(int(der_FARO), int(anc_FARO))>=int(sys.argv[5]):
                                                print chromosome+'\t'+coordinate+'\t'+var_ALPS+'\t'+var_FSCA+'\t'+var_FARO+'\t'+var_AMER
					else:
						if varlist_ALPS.index(min(varlist_ALPS)) == varlist_FSCA.index(min(varlist_FSCA)) == varlist_FARO.index(min(varlist_FARO)) == varlist_AMER.index(min(varlist_AMER)):
							pass
						else:
							print chromosome+'\t'+coordinate+'\t'+var_ALPS+'\t'+var_FSCA+'\t'+var_FARO+'\t'+var_AMER
		else:
			print 'ERROR - SITE DOES NOT CORRESPOND'


sys.exit()
