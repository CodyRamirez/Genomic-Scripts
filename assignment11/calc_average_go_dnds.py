#!/usr/bin/env python3
""" Finds all of the genes for each GO term in a GFF file and calculates the average dn/ds ratio for
that set of genes.

	Usage: 
		python3 calc_average_go_dnds.py <gff> <run_SNAP.py output file of dnds values>

		Example:
			python3 calc_average_go_dnds.py /home/assignments/assignment11/saccharomyces_cerevisiae.gff alignments_all_dnds.txt

	Output file name:
		average_go_dnds.txt
"""


import sys
import numpy as np
from collections import defaultdict


# sys.arg is a list containing 3 elements: the script name and 2 command line argument
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 3):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)

# Save the input arguments as varaibles
gff_file = open(sys.argv[1], 'r')
alignments_all_dnds_file = open(sys.argv[2], 'r')

gff_goid_gene_name_list = []
gene_name_list = []
dnds_ratio_dictionary = {}
gff_goid_gene_name_dictionary = defaultdict(list)



# Using a for loop to iterate through every line in alignments_all_dnds_file
for line in alignments_all_dnds_file:
	# Split line into list at whitespace, strip removes leading and trailing whitespace
	line_list = line.strip().split()

	# Using an IF statement to ensure that dn/ds header is skipped and all dh/ds ratios are converted to a float and added to the list
	if line_list[0] != 'Gene_names' and line_list[1] != 'dn/ds':

		# Appending all gene names to gene_name_list
		gene_name_list.append(line_list[0])

		# Appending all dn/ds ratio values to dnds_ratio_dictionary
		dnds_ratio_dictionary[line_list[0]] = float(line_list[1])


"""Parsing apart the GFF file to obtain a list of tuples containing the GO ID and its corresponding gene """
for line in gff_file:
	#
	line_list = line.strip().split()

	if not line.startswith("#") and line_list[2] == 'gene':
		
		for i in range(0, len(line_list[8].split(';'))):

			if 'Ontology_term' in line_list[8].split(';')[i]:

				for go_id in line_list[8].split(';')[i].split('=')[1].split(','):

					gene_name = line_list[8].split(';')[0].split('=')[1]

					if gene_name in gene_name_list:

						gff_goid_gene_name_list.append((go_id, gene_name))

"""Creating a dictionary where the key is the GO ID and the value is a list of all genes corresponding to the GO ID """
for go_id_tuple, gene_name_tuple in gff_goid_gene_name_list:

	gff_goid_gene_name_dictionary[go_id_tuple].append(gene_name_tuple)


# Opening the output file to write the output of the SNAP.pl program into
output_fileobject = open('average_go_dnds.txt', "w")
# Prints the header of the output file to the output file
print('GO ID:\tAverage dn/ds ratio:\tGene names:', file = output_fileobject)


"""Printing out columns corresponding to:
	0) GO ID
	1) Average dn/ds ratio of the genes annotated with the GO term
	2) Comma-separated list of gene names annotated with the GO term
"""
for key,value in gff_goid_gene_name_dictionary.items():

 	dnds_ratio_list = []

 	for gene_index in value:

 		dnds_ratio_list.append(dnds_ratio_dictionary[gene_index])

 	print(key + '\t' + str(np.mean(dnds_ratio_list)) + '\t' + str(value).strip('[').strip(']').replace("\'", ""), file = output_fileobject)

# Closing the output file
output_fileobject.close()