#!/usr/bin/env python3
"""
This script:
1) Takes a WGBS bed file and calculates the methylation level of each CpG
2) Writes a bed-like file of CpG methylation. The columns will be 1)chromosome, 2)start, 3)stop, and 4)methylation. Does not output CpGs that have 0X coverage. Saves the file as <WGBS bed basename>_CpG_methylation.bed
3) Plots the distribution of CpG methylation levels in a histogram as <WGBS bed basename>_methylation_distribution.png
4) Plots the distribution of read coverage for all CpGs for coverages between 0X and 100X as <WGBS bed basename>_CpG_coverage_distribution.png
5) Calculates and prints the fraction of CpGs that have 0X coverage

Usage: python3 analyze_WGBS_methylation.py <WGBS bed>
"""


# Import modules needed for program
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np


########## MAIN LOOP ##########

# sys.arg is a list containing 2 elements: the script name and 1 command line argument
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 2):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)

# Saves the input arguments as variables
BGM_WGBS_bed_file = open(sys.argv[1],"r")
path_to_file = sys.argv[1]

# Creating dynamic output filenames
# Grab just the filename
basename = os.path.basename(path_to_file)
# Remove the extension
basename_no_extension = os.path.splitext(basename)[0]

# Create output file name
output_file_name = basename_no_extension + "_CpG_methylation.bed"
# Open output file for writing
output_fileobject = open(output_file_name, "w")

# Initializing a CpG_methylation_level_list to store all CpG methylation values with great than zero coverage
CpG_methylation_level_list = []
# Initializing a CpG_coverage_list
CpG_coverage_list = []
# Initializing a CpG_zero_coverage_count
CpG_zero_coverage_count = 0
# Initializing a CpG_coverage_count
CpG_coverage_count = 0

# Use a for loop to iterate through every line in BGM_WGBS_bed_file
for line in BGM_WGBS_bed_file:
	line_list = line.strip().split()# Split line into list at whitespace, strip removes leading and trailing whitespace
	
	CpG_coverage = (float(line_list[3])+float(line_list[4]))# Calculating CpG_coverage
	CpG_coverage_list.append(CpG_coverage)# Appending the current CpG_coverage to CpG_coverage_list

	# This IF/ELSE is to ensure that NO CpGs with 0X coverage (0 # of C base calls and 0 # of T base calls) are output
	if line_list[3] == '0' and line_list[4] == '0':# Use a IF statement to determine if the number of C base calls and the number of T base calls equals zero then
		CpG_zero_coverage_count += 1# Increment CpG_zero_coverage_count by 1
	else:
		CpG_coverage_count += 1# Increment CpG_coverage_count by 1
		CpG_methylation_level = (float(line_list[3])/(float(line_list[3])+float(line_list[4])))# Calculates the methylation level for each CpG using the formula: CpG_methylation_level = (# C base calls)/(# C base calls + # T base calls)
		CpG_methylation_level_list.append(CpG_methylation_level)# Appends the current CpG_methylation_level value to the list
		# Prints chr, start, stop, and methylation level to an output file
		print(str(line_list[0]) + '\t' + str(line_list[1]) + '\t' + str(line_list[2]) + '\t' +  str(CpG_methylation_level), file = output_fileobject)

output_fileobject.close()# Closing the output file


# Create histogram plot name
histogram_plot_name = basename_no_extension + "_methylation_distribution.png"
fig = plt.figure()# Start a new figure
plt.hist(CpG_methylation_level_list, bins = 100, histtype = 'bar', align = 'mid')# PLotting CpG_methylation_level using 100 bins with bars aligned to the middle
plt.xlabel('CpG Methylation Level')# X-axis label
plt.ylabel('Counts')# Y-axis label
plt.title('Methylation Distribution in Chromosome 21')# Title of plot
fig.savefig(histogram_plot_name)# Saving the figure


# Create read coverage plot name
read_coverage_plot_name = basename_no_extension + "_CpG_coverage_distribution.png"
fig = plt.figure()# Start a new figure
plt.hist(CpG_coverage_list, bins = 100, range = (0,100), histtype = 'bar', align = 'mid')# Plotting CpG_coverage using 100 bins with bars aligned to the middle
plt.xlabel('Coverage')# X-axis label
plt.ylabel('Counts')# Y-axis label
plt.title('Coverage Distribution in Chromosome 21')# Title of plot
fig.savefig(read_coverage_plot_name)# Saving the figure


# Prints the fraction of CpGs that have zero coverage
print('The fraction of CpGs that have 0X coverage is: ' + str(CpG_zero_coverage_count/(CpG_zero_coverage_count+CpG_coverage_count)))
