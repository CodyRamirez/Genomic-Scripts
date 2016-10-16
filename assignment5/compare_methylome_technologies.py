#!/usr/bin/env python3
"""
This script:
1) Creates the following scatter plots:
	a) MeDIP-seq RPKM vs. MRE-seq RPKM. Saved the plot as <MeDIP-seq RPKM bed basename>_vs_<MRE-seq RPKM bed basename>.png
	b) MeDIP-seq RPKM vs. WGBS average DNA methylation level. Save the plot as <MeDIP-seq RPKM bed basename>_vs_<WGBS methylation level bed basename>.png
	c) MRE-seq RPKM vs. WGBS average DNA methylation level. Save the plot as <MRE-seq RPKM bed basename>_vs_<WGBS methylation level bed basename>.png
2) Calculates the correlation for each comparision.
3) Prints the correlation to stdout

Usage: python3 compare_methylome_technologies.py <MeDIP-seq RPKM bed> <MRE-seq RPKM bed> <WGBS methylation level bed>
"""


# Import modules needed for program
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


########## MAIN LOOP ##########

# sys.arg is a list containing 4 elements: the script name and 1 command line argument
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 4):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)

# Saves the input arguments as variables
MeDIP_CGI_RPKM_file = open(sys.argv[1],"r")
MRE_CGI_RPKM_file = open(sys.argv[2], "r")
WGBS_CGI_methylation_file = open(sys.argv[3], "r")
path_to_MeDIP_CGI_RPKM_file = sys.argv[1]
path_to_MRE_CGI_RPKM_file = sys.argv[2]
path_to_WGBS_CGI_methylation_file = sys.argv[3]

# Creating dynamic output filenames
# Grab just the filename
MeDIP_basename = os.path.basename(path_to_MeDIP_CGI_RPKM_file)
MRE_basename = os.path.basename(path_to_MRE_CGI_RPKM_file)
WGBS_basename = os.path.basename(path_to_WGBS_CGI_methylation_file)

# Remove the extension
MeDIP_basename_no_extension = os.path.splitext(MeDIP_basename)[0]
MRE_basename_no_extension = os.path.splitext(MRE_basename)[0]
WGBS_basename_no_extension = os.path.splitext(WGBS_basename)[0]

# Initializing list to later contain the values of the files
MeDIP_CGI_RPKM_list = []
MRE_CGI_RPKM_list = []
WGBS_CGI_methylation_list = []


# Use a for loop to iterate through every line in CGI_methylation_file
for line in MeDIP_CGI_RPKM_file:
	line_list = line.strip().split()# Split line into list at whitespace, strip removes leading and trailing whitespace
	if line_list[1] not in ('15035943', '15038195', '15052974', '15068825', '15077051', '9825442'):# Excluding lines that are not shared by WGBS_CGI_methylation_file and the outlier
		MeDIP_CGI_RPKM_list.append(float(line_list[4]))# Appending all values of MeDIP_CGI_RPKM that pass the filter to the list
# Use a for loop to iterate through every line in MRE_CGI_methylation_file
for line in MRE_CGI_RPKM_file:
	line_list = line.strip().split()# Split line into list at whitespace, strip removes leading and trail whitespace
	if line_list[1] not in ('15035943', '15038195', '15052974', '15068825', '15077051', '9825442'):# Excluding lines that are not shared by WGBS_CGI_methylation_file and the outlier
		MRE_CGI_RPKM_list.append(float(line_list[4]))# Appending all values of MRE_CGI_RPKM_file that pass the filter to the list
# Use a for loop to iterate through every line in WGBS_CGI_methylation_file
for line in WGBS_CGI_methylation_file:
	line_list = line.strip().split()# Split line into list at whitespace, strip removes leading and trailing whitespace
	if line_list[1] != '9825442':# Excluding the outlier
		WGBS_CGI_methylation_list.append(float(line_list[3]))# Appending all values of WGBS_CGI_methylation that pass the filter to the list


first_scatter_plot_name = MeDIP_basename_no_extension + "_vs_" + MRE_basename_no_extension + "_outliers_removed.png"# Create first scatter plot name
fig = plt.figure()# Start a new figure
plt.scatter(MeDIP_CGI_RPKM_list, MRE_CGI_RPKM_list)# Plotting MeDIP_CGI_RPKM and MRE_CGI_RPKM against each other using a scatter plot
plt.xlabel('MeDIP-seq RPKM')# X-axis label
plt.ylabel('MRE-seq RPKM')# Y-axis label
plt.title('MeDIP-seq RPKM vs. MRE-seq RPKM')# Title of scatter plot
fig.savefig(first_scatter_plot_name)# Saving the figure
print('The correlation between MeDIP-seq RPKM and MRE-seq RPKM is:' + str(stats.spearmanr(MeDIP_CGI_RPKM_list, MRE_CGI_RPKM_list)))# Printing out the Spearman's correlation for MeDIP_CGI_RPLM vs. MRE_CGI_RPKM



second_scatter_plot_name = MeDIP_basename_no_extension + "_vs_" + WGBS_basename_no_extension + "_outliers_removed.png"# Create second scatter plot name
fig = plt.figure()# Start a new figure
plt.scatter(MeDIP_CGI_RPKM_list, WGBS_CGI_methylation_list)# Plotting MeDIP_CGI_RPKM and WGBS_CGI_methylation against each other using a scatter plot
plt.xlabel('MeDIP-seq RPKM')# X-axis label
plt.ylabel('WGBS methylation')# Y-axis label
plt.title('MeDIP-seq RPKM vs. WGBS methylation')# Title of scatter plot
fig.savefig(second_scatter_plot_name)#Saving figure
print('The correlation between MeDIP-seq RPKM and WGBS methylation is:' + str(stats.spearmanr(MeDIP_CGI_RPKM_list, WGBS_CGI_methylation_list)))# Printing out the Spearman's correlation for MeDIP_CGI_RPKM vs. WGBS_CGI_methylation


thrid_scatter_plot_name = MRE_basename_no_extension + "_vs_" + WGBS_basename_no_extension + "_outliers_removed.png"# Create third scatter plot name
fig = plt.figure()# Start a new figure
plt.scatter(MRE_CGI_RPKM_list, WGBS_CGI_methylation_list)# Plotting MRE_CGI_RPKM and WGBS_CGI_methylation against each other using a scatter plot
plt.xlabel('MRE-seq RPKM')# X-axis label
plt.ylabel('WBGS methylation')# Y-axis label
plt.title('MRE-seq RPKM vs. WBGS methylation')# Title of scatter plot
fig.savefig(thrid_scatter_plot_name)# Saving figure
print('The correlation between MRE-seq RPKM and WGBS methylation is:' + str(stats.spearmanr(MRE_CGI_RPKM_list, WGBS_CGI_methylation_list)))# Print out the Spearman's correlation for MRE_CGI_RPKM vs. WGBS_CGI_methylation
