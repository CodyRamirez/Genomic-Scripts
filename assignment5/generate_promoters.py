#!/usr/bin/env python3
"""
This script:
1) Takes a bed file of gene coorinates and creates a bed file of their promoters

Definition of promoters within this code:
1) On average promoters are approximately 250 base pairs upstream of the start site
2) Promoters can vary from 100 - 1,000 base pairs, therefore I am using 500 base pairs as my promoter lengths

Output: The columns will be 1) chromosome, 2) promoter start position, 3) promoter stop position, 4) gene name, 5) strand
Usage: python3 generate_promoters.py <bed of gene coordinates>
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
gene_coordinates_file = open(sys.argv[1],"r")
path_to_file = sys.argv[1]

# Creating dynamic output filenames
# Grab just the filename
basename = os.path.basename(path_to_file)
# Remove the extension
basename_no_extension = os.path.splitext(basename)[0]

# Create output file name
output_file_name = basename_no_extension + "_promoters.bed"
# Open output file for writing
output_fileobject = open(output_file_name, "w")


# Use a for loop to iterate through every line in gene_coordinate_file
for line in gene_coordinates_file:
	line_list = line.strip().split()# Split line into list at whitespace, strip removes leading and trailing whitespace
	promoter_start = int(line_list[1]) - 750# Creating the start position of the promoter by substracting 750 from the start position of the refGene
	promoter_stop = int(line_list[1]) - 250# Creating the stop position of the promoter by substracting 250 from the start position of the refGene
	print(str(line_list[0]) + '\t' + str(promoter_start) + '\t' + str(promoter_stop) + '\t' +  str(line_list[3]) + '\t' + str(line_list[5]), file = output_fileobject)# Print to the output file
output_fileobject.close()# Closing the output file
