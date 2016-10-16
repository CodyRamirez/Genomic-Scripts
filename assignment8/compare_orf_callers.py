#!/usr/bin/env python3
"""Compares the predicted ORFs from call_orfs.py script to those called by MetaGeneMark.

Usage: python3 compare_orf_callers.py <fasta 1> <fasta 2>

Args:
	fasta 1 = Fasta file outputs from call_orfs.py
	fasta 2 = Fasta file outputs from MetaGeneMark

Outputs:
	1) ORFs found in both the outputs of call_orfs.py and MetaGeneMark
	2) ORFs unique to the output of call_orfs.py
	3) ORFs unique to the output of MetaGeneMark
	4) Number of ORFs in each of these three categories
"""


import sys

# sys.arg is a list containing 3 elements: the script name and 2 command line argument
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 3):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)


# Save the input arguments as varaibles
call_orfs_fasta_file = open(sys.argv[1], "r")
meta_gene_mark_fasta_file = open(sys.argv[2], "r")

call_orfs_dictionary = {}
meta_gene_mark_dictionary = {}

# Use a for loop to iterate through the file line by line
for line in call_orfs_fasta_file:
	# Remove newline
	line = line.rstrip()
	# (Getting the sequence name) If the line starts with the character ">" then,
	if line.startswith(">"):
		# Assign seq_name
		seq_name = line
	# (Read in sequence in uppercase) If the line does not start with the character ">" then,
	if not line.startswith(">"):
		# Make all of the characters within the line uppercase
		line = line.upper()
	# Put the seq_name and sequence in call_orfs_dictionary
	call_orfs_dictionary[seq_name] = line



# Use a for loop to iterate through the file line by line
for line in meta_gene_mark_fasta_file:
	# Remove newline
	line = line.strip()
	seq_line = ''
	# (Getting the sequence name) If the line starts with the character ">" then,
	if line.startswith(">"):
		# Assign seq_name
		seq_name = line
	# (Read in sequence in uppercase) If the line does not start with the character ">" then,
	if not line.startswith(">"):
		# Make all of the characters within the line uppercase
		line = line.upper()
	# Put the seq_name and sequence in meta_gene_mark_dictionary
	meta_gene_mark_dictionary[seq_name] = line



found_in_both_dictionary = {}
found_in_call_orfs_dictionary = {}
found_in_meta_gene_mark_dictionary = {}

for call_orfs_name, call_orfs_sequence in call_orfs_dictionary.items():
	for meta_gene_mark_name, meta_gene_mark_sequence in meta_gene_mark_dictionary.items():
		if call_orfs_sequence == meta_gene_mark_sequence:
			found_in_both_dictionary[call_orfs_name+"_"+meta_gene_mark_name] = call_orfs_sequence


for call_orfs_name, call_orfs_sequence in call_orfs_dictionary.items():
	for meta_gene_mark_name, meta_gene_mark_sequence in meta_gene_mark_dictionary.items():
		if call_orfs_sequence != meta_gene_mark_sequence and call_orfs_sequence not in found_in_both_dictionary.values() and meta_gene_mark_sequence not in found_in_both_dictionary.values():
			found_in_call_orfs_dictionary[call_orfs_name] = call_orfs_sequence
			found_in_meta_gene_mark_dictionary[meta_gene_mark_name] = meta_gene_mark_sequence



print('\nORFs found in both')
for found_in_both_name, found_in_both_sequence in found_in_both_dictionary.items():
	print(found_in_both_name+'\n'+found_in_both_sequence)

print('\nORFs found only in call_orfs.py')
for found_in_call_orfs_name, found_in_call_orfs_sequence in found_in_call_orfs_dictionary.items():
	print(found_in_call_orfs_name+'\n'+found_in_call_orfs_sequence)

print('\nORFs found only in MetaGeneMark')
for found_in_meta_gene_mark_name, found_in_meta_gene_mark_sequence in found_in_meta_gene_mark_dictionary.items():
	print(found_in_meta_gene_mark_name+'\n'+found_in_meta_gene_mark_sequence)

print('The number of ORFs found in both: '+str(len(found_in_both_dictionary)))
print('The number of ORFs found only in call_orfs.py: '+str(len(found_in_call_orfs_dictionary)))
print('The number of ORFs found only in MetaGeneMark: '+str(len(found_in_meta_gene_mark_dictionary)))

