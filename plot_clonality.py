#!/usr/bin/env python3
"""
Plots the clonality of two samples based on VAF.

Usage: python3 plot_clonality.py <bam_readcount> <xlabel> <xVAF> <ylabel> <yVAF>

Required Arguments:
    <bam_readcount> = Path to bam-readcount output file
    <xlabel> = X-axis label
    <xVAF> = Column number containing the VAF numbers for x-axis data
    <ylabel> = Y-axis label
    <yVAF> = Column number containing the VAF numbers for y-axis data
"""

# Import modules needed for program
import sys
import matplotlib.pyplot as plt

# sys.arg is a list containing 6 elements: 0) script name, 1) bam_readcount_output, 2) xlabel, 3) xVAF, 4) ylabel, 5) yVAF
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 6):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)

# Saves the input arguments as variables
bam_readcount_output_file = open(sys.argv[1],"r")
x_label = sys.argv[2]
x_VAF_column_number = int(sys.argv[3]) - 1
y_label = sys.argv[4]
y_VAF_column_number = int(sys.argv[5]) - 1

# Initializing list to contain the VAFs
x_VAF_numbers = []
y_VAF_numbers = []

# Skipping the first line in the file to ensure header is removed
next(bam_readcount_output_file)

# Use a for loop to iterate through every line in the bam-readcount output file
for line in bam_readcount_output_file:
	line = line.strip().split()# Split line into list at whitespace, strip removes leading and trailing whitespace
	if line[x_VAF_column_number] != 'NA' and line[y_VAF_column_number] != 'NA':
		x_VAF_numbers.append(float(line[x_VAF_column_number]))
		y_VAF_numbers.append(float(line[y_VAF_column_number]))
	else:
		print(line)

figure_name = x_label + '_' + y_label + '_clonality_plot.png'
figure_title = x_label + ' vs. ' + y_label
fig = plt.figure()# Start a new figure
plt.scatter(x_VAF_numbers, y_VAF_numbers)#
plt.xlim(0,100)
plt.ylim(0,100)
plt.xlabel(x_label)# X-axis label
plt.ylabel(y_label)# Y-axis label
plt.title(figure_title)# Title of scatter plot
fig.savefig(figure_name)# Saving the figure