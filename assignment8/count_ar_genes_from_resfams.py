#!/usr/bin/env python3
"""Takes resfams_annotations.txt as input and calculates the number of antibiotic resistance genes identified
	in contigs

Usage: python3 count_ar_genes_from_resfams.py <hmmscan output>

Args:
	hmmscan output = Path to hmmscan output

Outputs:
	The number of antibiotic resistance genes identified in contigs
"""

import sys

# sys.arg is a list containing 2 elements: the script name and 1 command line argument
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 2):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)

# Save the input argument as contigs_file varaible
resfams_annotations_file = open(sys.argv[1], "r")
# Intializing a variable count
count = 0
# Use a for loop to iterate through the file line by line
for line in resfams_annotations_file:
	# Remove newline
	line = line.rstrip()
	# Keeping track of the number of genes
	count += line.count('gene_id')
# Printing out the number of genes
print('The number of antibiotic resistance genes identified in contigs: '+str(count))