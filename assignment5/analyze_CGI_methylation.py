#!/usr/bin/env python3
"""
This script:
1) Takes the average CGI methylation bed file and plots the distribution of average CGI methylation levels.

Usage: python3 analyze_CGI_methylation.py <average CGI methylation bed>
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
WGBS_CGI_methylation_file = open(sys.argv[1],"r")
path_to_file = sys.argv[1]

# Creating dynamic output filenames
# Grab just the filename
basename = os.path.basename(path_to_file)
# Remove the extension
basename_no_extension = os.path.splitext(basename)[0]

# Initializing average_WGBS_CGI_methylation_list
average_WGBS_CGI_methylation_list = []


# Use a for loop to iterate through every line in WGBS_ CGI_methylation_file
for line in WGBS_CGI_methylation_file:
	line_list = line.strip().split()# Split line into list at whitespace, strip removes leading and trailing whitespace
	average_WGBS_CGI_methylation_list.append(float(line_list[3]))# Appending the current average WGBS CGI methylation value to average_WGBS_CGI_methylation_list


histogram_plot_name = basename_no_extension + ".png"# Create histogram plot name
fig = plt.figure()# Start a new figure
plt.hist(average_WGBS_CGI_methylation_list, bins = 100, histtype = 'bar', align = 'mid')# Plotting average WGBS CGI methylation using 100 bins with bars aligned to the middle
plt.xlabel('Average CGI Methylation Level')# X-axis label
plt.ylabel('Counts')# Y-axis label
plt.title('Average Promoter CGI Methylation Distribution in Chromosome 21')# Title of plot
fig.savefig(histogram_plot_name)# Saving figure
