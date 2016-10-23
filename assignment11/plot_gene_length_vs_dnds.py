#!/usr/bin/env python3
"""	Creates a scatterplot of the gene lengths (y-axis) versus the dn/ds ratio (x-axis):
		1) Restrict the range of the x-axis to avoid distortion caused by outliers
		2) Save the plot as gene_length_vs_dnds.png
		3) Add informative axis labels, a title, etc.

	Usage: 
		python3 plot_gene_length_vs_dnds.py <run_SNAP.py error file> <run_SNAP.py output file of dnds values>

		Example:
			python3 plot_gene_length_vs_dnds.py alignments.err alignments_all_dnds.txt
"""


import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



# sys.arg is a list containing 3 elements: the script name and 2 command line argument
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 3):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)

# Save the input arguments as varaibles
alignments_err_file = open(sys.argv[1], 'r')
alignments_all_dnds_file = open(sys.argv[2], 'r')

# Initializing list to hold all gene length values
gene_length_list = []
# Initializing list to hold all dn/ds ratio values
dnds_ratio_list = []
# Initializing a dictionary to contain gene names with it's corresponding gene lengths
gene_length_dictionary = {}


# Using a for loop to iterate through every line in alignments_err_file
for line in alignments_err_file:
	# Split line into list at whitespace, strip removes leading and trailing whitespace
	line_list = line.strip().split()

	# Filtering out and adding gene names with their corresponding gene lengths
	if line_list[1] != 'ERROR:':
		if line_list[0] == 'Sequence':
			gene_length = line_list[2]
			gene_length_dictionary[gene_name] = gene_length
		else:
			gene_name = line_list[0]


# Using a for loop to iterate through every line in alignments_all_dnds_file
for line in alignments_all_dnds_file:
	# Split line into list at whitespace, strip removes leading and trailing whitespace
	line_list = line.strip().split()

	if line_list[0] != 'Gene_names':
		# Appending all gene length values to gene_length_list
		gene_length_list.append(int(gene_length_dictionary[line_list[0]]))

	# Using an IF statement to ensure that dn/ds header is skipped and all dh/ds ratios are converted to a float and added to the list
	if line_list[1] != 'dn/ds':
		# Appending all dn/ds ratio values to dnds_ratio_list
		dnds_ratio_list.append(float(line_list[1]))



# Starting a new figure
fig = plt.figure()
# Restricting the range of the x-axis to avoid distortion caused by outliers
plt.xlim(0,0.6)
# Setting y-axis to avoid default spacing 
plt.ylim(0,1400)
# Plotting gene lengths and dn/ds ratios against each other in a scatter plot
plt.scatter(dnds_ratio_list, gene_length_list)
# Labeling the X-axis
plt.xlabel('dn/ds Ratio Values')
# Labeling the Y-axis
plt.ylabel('Gene Length Values')
# Initializing the title of the scatter plot
plt.title('Gene Lengths vs. dn/ds Ratios')
# Saving the figure as gene_length_vs_dnds.png
fig.savefig('gene_length_vs_dnds.png')