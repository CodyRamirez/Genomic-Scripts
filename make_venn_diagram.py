#!/usr/bin/env python3
"""
Creates a venn diagram from a tab-delimited file. Using specified columns as unique identifiers and 
using other specified columns to count the number of occurrences.

    Usage: 
    	python3 make_venn_diagram.py <file> <uic> <vc> <ds1> <ds_vc1> <ds2> <ds_vc2> <ds3> <ds_vc3>

    Example:
    	python3 make_venn_diagram.py bam_readcount_H_VQ_GTB22.tsv 1-5 C001 11 L001 14 L002 17
    OR
    	python3 make_venn_diagram.py bam_readcount_H_VQ_GTB22.tsv 1:5 C001 11 L001 14 L002 17
    
    Required Arguments:
        <file> = Full path to file
        <uic> = Requires TWO numbers to mark the beginning and end of unique identifying columns to be used
        <vc> = Minimum variant count required in order to add to data set list
        <ds1> = Label for data set 1
        <ds_vc1> = Column number containing the variant count numbers for data set 1
        <ds2> = Label for data set 2
        <ds_vc2> = Column number containing the variant count numbers for data set 2
    
    Optional Arguments:
        <ds3> = Label for data set 3
        <ds_vc3> = Column number containing the variant count numbers for data set 3

    Note:
    	I assume a single line header exist and skip it (line 57)
"""

# Import modules needed for program
import sys, matplotlib, matplotlib_venn

# sys.arg is a list containing at least 6 elements: 0) script name, 1) path to file, 2) data_set_1_label, 3) data_set_1, 4) data_set_2_label, 5) data_set_2; Optional 6) data_set_3_label, 7) data_set_3
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) <= 6):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)

# Saving all inputs as variable to use later
read_file = open(sys.argv[1],"r")
beginning = int(sys.argv[2][0]) - 1
end = int(sys.argv[2][-1])
var_count = int(sys.argv[3])
data_set_1_label = sys.argv[4]
data_set_1 = int(sys.argv[5]) - 1
data_set_2_label = sys.argv[6]
data_set_2 = int(sys.argv[7]) - 1
if (len(sys.argv) == 10):
    data_set_3_label = sys.argv[8]
    data_set_3 = int(sys.argv[9]) - 1

# Initializing list to contain a list of all unique identifiers present within a data set
ds_1_var_count = set()
ds_2_var_count = set()
ds_3_var_count = set()
# Skipping the first line, assumes a header exist
next(read_file)
# Reading in the file a line at a time
for line in read_file:
	line = line.split()
	# Saving the list of unique identifiers as a single string
	unique_var = '_'.join(line[beginning:end])
	# Removing any variants with NA values
	if 'NA' not in line:
		# Only adding unique identifiers with a variant count greater or equal to a set value for all data sets
		if int(line[data_set_1]) >= var_count:
			ds_1_var_count.add(unique_var) 
		if int(line[data_set_2]) >= var_count:
			ds_2_var_count.add(unique_var)
		if (len(sys.argv) == 10):
			if int(line[data_set_3]) >= var_count:
				ds_3_var_count.add(unique_var)
# Plotting a venn diagram with the given information
fig = matplotlib.pyplot.figure()
if (len(sys.argv) == 10):
	matplotlib_venn.venn3([ds_1_var_count, ds_2_var_count, ds_3_var_count], set_labels = (data_set_1_label, data_set_2_label, data_set_3_label))
	fig.savefig(data_set_1_label + '_' + data_set_2_label + '_' + data_set_3_label + '_venn_diagram.pdf', format = 'pdf')
else:
	matplotlib_venn.venn2([ds_1_var_count, ds_2_var_count], set_labels = (data_set_1_label, data_set_2_label))
	fig.savefig(data_set_1_label + '_' + data_set_2_label + '_venn_diagram.png', format = 'pdf')
