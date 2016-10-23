#!/usr/bin/env python3
"""	1) Runs SNAP.pl on every fasta file in a given directory
	2) Prints the substitution rates between Scer (S. cerevisiae) and Spar (S. paradoxus)
from each fasta file to an output file.
		a) DO NOT include genes with dn/ds values of 0 or NA in your output
		b) Output file format: Tab-delimited columns will be gene name and dn/ds ratio.
		c) Save the file as <alignments directory basename>_all_dnds.txt
			i) Automatically parse out the basename of the alignments directory
	3) Prints the average dn/ds ratio to standard output
	4) Prints the genes and their dn/ds ratio for genes with dn/ds > 1 to standard output

	Usage: 
		python3 run_SNAP.py <alignments directory> <output directory> 2> <error file>

		Example:
			python3 run_SNAP.py /home/assignments/assignment11/alignments part2/ 2> alignments.err 
"""


import sys, os, subprocess
import numpy as np


# sys.arg is a list containing 3 elements: the script name and 2 command line argument
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 3):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)

# Save the input arguments as varaibles
alignment_directory = sys.argv[1]
output_directory = sys.argv[2]


# Initializing a list containing the names of the enteries in the directory given by path
file_name = os.listdir(alignment_directory)
# Initializing a list to hold the full path name of the files for easy access
full_path_to_file_list = []
# Iterating through the list of file_names to create a list of full path names to the file
for i in file_name:
	# Appending full path names of the files to a list using os.path.join to create a full path name given a directory and the file names within it
	full_path_to_file_list.append(os.path.join(alignment_directory, i))


# Creating dynamic output filenames
# Grabbing the basename of the alignment directory
basename = os.path.basename(alignment_directory)
# Initializing output_file_name to the alignment directory basename and a descriptor
output_file_name = basename + "_all_dnds.txt"
# Opening the output file to write the output of the SNAP.pl program into
output_fileobject = open(output_file_name, "w")
# Prints the header of the output file to the output file
print('Gene_names\t dn/ds', file = output_fileobject)


# Iterating through the list of full_path_to_file_list to obtain the fasta_file_path
for fasta_file_path in full_path_to_file_list:
	# Writing the gene name to alignment.err file
	sys.stderr.write(os.path.basename(fasta_file_path).split('.')[0] + '\t')
	# Flushing the buffer after writing to alignment.err to ensure that each gene is written to it's corresponding run
	sys.stderr.flush()
	# Creating a list of commands to run within subprocess
	cmd = ['perl', '/home/assignments/assignment11/SNAP.pl', fasta_file_path, output_directory]
	# Running run_SNAP.pl on every fasta file within the directory /home/assignments/assignment11/alignments and outputting the resulting files into the designated output_directory
	subprocess.call(cmd)

	


# Initializing a list containing the names of the enteries in the output directory: part2/
alignment_file_names = os.listdir(output_directory)
# Initialzing a list to hold all dn/ds values found within output directory: part2/ 
dnds_value_list = []

# Iterating through current_file_name which is the current file name within the output directory of the fasta alignment: part2/
for current_file_name in alignment_file_names:
	# Initializing alignment_file to the the full path to the current file within the output directory: part2/
	alignment_file = os.path.join(output_directory, current_file_name)
	# Opening the current alignment_file output to read through the lines of the file. This would allow access to dn/ds ratios
	open_alignment_file = open(alignment_file, 'r')

	# Using a for loop to iterate through all lines within the current alignment_file output
	for line in open_alignment_file:
		# Parsing the current line within the file, stripping out all whitespace, then splitting on the whitespace and initializing all values to line_list
		line_list = line.strip().split()

		# Using an IF statement to ensure that the current line contains 'Scer' and 'Spar' but dn/ds ratio is greater than 0 and exist
		if 'Scer' in line_list and 'Spar' in line_list and line_list[12] != 'NA' and float(line_list[12]) > 0:
			# Spliting the current_file_name by a '.' and initializing the first element to gene_name
			gene_name = current_file_name.split('.')[0]
			# Printing out the gene_name and it's corresponding dn/ds ratio value alignment between 'Scer' and 'Spar' to alignments_all_dnds.txt
			print(gene_name + '\t' + line_list[12], file = output_fileobject)
			
			# Appending all dn/ds ratios value alignments between 'Scer' and 'Spar' to dnds_value_list
			dnds_value_list.append(float(line_list[12]))
			# Using an IF statement to filter for dn/ds ratio values greater than one then,
			if float(line_list[12]) > 1:
				# Printing out the gene_name and it's corresponding dn/ds ratio value to standard output
				print(gene_name + '\t' + line_list[12])


# Initializing the average dn/ds ratio to the variable average_dnds_value
average_dnds_value = np.mean(dnds_value_list)
# Printing out the average dn/ds ratio value to standard output
print('The average dn/ds ratio is: ' + str(average_dnds_value))



# Closing the output file
output_fileobject.close()