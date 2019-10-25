# J. Melo-Ferreira 10/04/2018
##
# Combines chr|site keys from different populations in a non-duplicated and ordered list
##
# Usage:
# python 1-combine-keys.py keys-ALPS keys-FSCA keys-FARO keys-AMER > output
##

import sys

# ALPS

ALPS_key_file = open(sys.argv[1], 'r')
ALPS_key_list = ALPS_key_file.readlines()

ALPS_key_file.close()



# ADD FSCA

FSCA_key_file = open(sys.argv[2], 'r')
FSCA_key_list = FSCA_key_file.readlines()

FSCA_key_file.close()



# COMBINE ALPS + FSCA

combined_list = ALPS_key_list + FSCA_key_list

ALPS_key_list = []
FSCA_key_list = []


combined_list = set(combined_list)
combined_list = list(combined_list)


# ADD FARO

FARO_key_file = open(sys.argv[3], 'r')
FARO_key_list = FARO_key_file.readlines()

FARO_key_file.close()

combined_list = combined_list + FARO_key_list

FARO_key_list = []

combined_list = set(combined_list)
combined_list = list(combined_list)



# ADD AMER

AMER_key_file = open(sys.argv[4], 'r')
AMER_key_list = AMER_key_file.readlines()

AMER_key_file.close()

combined_list = combined_list + AMER_key_list

AMER_key_list = []

combined_list = set(combined_list)
combined_list = list(combined_list)



# ORDER KEYS

combined_list = sorted(combined_list, key = lambda x: (x.split('|')[0], int(x.split('|')[1])))



print ''.join(combined_list)



sys.exit()
