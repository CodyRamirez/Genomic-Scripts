#!/usr/bin/env python3
"""Catenates "ref_count var_count VAF" columns of all files in a directory with .tsv into a single file as below:
	"chr start stop ref var ref_count_1 var_count_1 VAF_1 ... ref_count_N var_count_N VAF_N"

	Usage: 
		python3 pon_readcount_merge_v2.py <data directory> <existing directory path and/or output file name>

	Example:
		python3 pon_readcount_merge_v2.py lyphoma_data processed_data/all_variants_merged.tsv
"""
import sys, os, glob, pandas, numpy
# sys.arg is a list containing 3 elements: The script, directory where the data is stored, desired name for output file
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 3):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)
# Save the input arguments as varaibles
data_directory = sys.argv[1]
output_file_name = sys.argv[2]
# Making sure the user is not outputting to the same location as the data data being processed
if '/' in output_file_name and data_directory == output_file_name.split('/')[-2]:
	sys.exit("ERROR: You can not save the output in the same directory being processed.\nThis will cause the program to read the empty output file and crash.\nPlease pick a new data directory for output, thank you.")
# Setting a boolean statement to read in a single file
first_file = True
# Using a foor loop to interate through all files ending with '.tsv' within a given directory
for filename in glob.glob(os.path.join(data_directory, '*.tsv')):
	# Reading in the current data file as a DataFrame
	open_file_df = pandas.read_table(filename)
	# Grabbing the base name of the file to maintain column identity
	base_name = filename.split('/')[1].split('r')[0]
	# Renaming ref_count, var_count and VAF column values to maintain data identity 
	for x in range(5,8):
		open_file_df.rename(columns = {list(open_file_df)[x] : base_name+list(open_file_df)[x]}, inplace = True)
	# Removing all duplicates
	open_file_df.drop_duplicates(inplace = True)
	# Identifing the first line
	if first_file == True:
		first_file = False
		# Using the first file read in to seed the output
		seed_file = open_file_df
	else:
		# Merging all other files to the seed_file based on comparison to variant info columns 0-4 (chromosome_name, start, stop, reference, varinat)
		seed_file = pandas.merge(seed_file, open_file_df, how = 'left', on = list(seed_file)[0:5])
		# Checking for any null values in the DataFrame caused by variant rows not matching
		if seed_file.isnull().values.any() == True:
			sys.exit("ERROR: File " + base_name + " has a merge error!\nRow " + str(numpy.where(pandas.isnull(seed_file))[0][0]) + " does not match all variant info fields.\n\nSeed file variant info:\n" + str(seed_file.ix[numpy.where(pandas.isnull(seed_file))[0][0],0:5]) + "\n\n" + base_name + " file variant info:\n" + str(open_file_df.ix[numpy.where(pandas.isnull(seed_file))[0][0], 0:5]) + "\n\nPlease address this issue before proceeding!")
# Outputting the DataFrame as a tab-delimited file
seed_file.to_csv(path_or_buf = output_file_name, sep = '\t', index = False)