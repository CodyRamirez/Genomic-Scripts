#!/usr/bin/env python3
""" Catenates "ref_count var_count VAF" columns of all files in a directory with .tsv into a single file as below:
	"chr start stop ref var ref_count_1 var_count_1 VAF_1 ... ref_count_N var_count_N VAF_N"
	Usage: 
		python3 pon_readcount_merge.py <data directory> <output file name>
		Example:
			python3 pon_readcount_merge.py lyphoma_data all_variants_merged.tsv
"""
import sys, os, glob
# sys.arg is a list containing 3 elements: The script, directory where the data is stored, desired name for output file
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 3):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)
# Save the input arguments as varaibles
data_directory = sys.argv[1]
output_file_name = sys.argv[2]
# Opening an output file to write the merged output to
output_fileobject = open(output_file_name, "w")
# Setting a boolean statement to read in a single file
first_file = True
# Initializing list to contain the header and variant output
header = []
variant_list = []
# Using a foor loop to interate through all files ending with '.tsv' within a given directory
for filename in glob.glob(os.path.join(data_directory, '*.tsv')):
	# Initializing line_number variable to use downstream for exact row comparison
	line_number = 0
	# Opening all files as read only to access data
	open_file = open(filename, "r")
	# Sorting and deduplicating files
	open_file = sorted(list(set(open_file)))
	# Grabbing the base name of the file to maintain column identity
	base_name = filename.split('/')[1].split('r')[0]
	# Only allowing one file in to seed the output file
	if first_file == True:
		first_file = False
		# Iterating through every line in the open files
		for line in open_file:
			line = line.split()
			# Checking for the header
			if line[0][0].isalpha() == True:
				for x in range(5,8):
					# Adding the file names to the data columns in the header 
					line[x] = base_name+line[x]
				# Saving formatted header
				header = line
			else:
				# Appending all other variant lines to the variant_list
				variant_list.append(line)
	# Processing the rest of the files
	else:
		# Iterating through every line in the open files
		for line in open_file:
			line = line.split()
			# Checking for the header
			if line[0][0].isalpha() == True:
				for y in range(5,8):
					# Adding the file names to the data columns in the header
					header.append(base_name+line[y])
			else:
				# Comparing the current row to the exact corresponding row number within the files to maintain specificity and avoid another for loop
				if line[0:5] == variant_list[line_number][0:5]:
					for z in range(5,8):
						# Appends the new values 'ref_count var_count VAF' for a given file 
						variant_list[line_number].append(line[z])
					# Incrementing line_number to maintain current with the row index
					line_number += 1
				else:
					sys.exit("ERROR: Line " + str(line_number) + " within the seed file does not match the variant info columns.\nSeed file columns: " + str(variant_list[line_number][0:5]) + "\nCurrent file (" + filename + ") columns: " + str(line[0:5]) + "\nPlease address discrepancy and restart program.")

# Prints header and corresponding values for variants
print(*header, sep = '\t', file = output_fileobject)
for variant_out in variant_list:
	print(*variant_out, sep = '\t', file = output_fileobject)
# Closing the output file
output_fileobject.close()
