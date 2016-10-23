#!/usr/bin/env python3
"""	1) Calculates the mean in variant coverage for each dataset
	2) Calculates the variance in variant coverage for each dataset
	3) Calculates the number of variants with unknown coverages for each dataset
	4) In a SINGLE figure, makes a boxplot for each dataset to show the coverage distribution
	5) Using the heterozygous sites identified in the blood_standard dataset,
		Calculates the rate of allelic dropout for each of the remaining datasets
			a) Report the numerator and denominator used to calculate this rate
	6) Finds the set of "live" variants AND "dead" variants
		I) Calculates and prints the number of variants in each set to STDOUT
		II) Calculates and prints to STDOUT the number of:
			a) Transitions for each 
			b) Transversions
			c) Transition/Transversion ratio


	Usage: 
		python3 genotype_analysis.py --vcf <VCF>

		Example:
			python3 genotype_analysis.py --vcf sperm_genotype_calls.vcf

	Args:
		VCF = The VCF file, sperm_genotype_calls.vcf, contains genotype calls and additional
data summaries for 9 sequencing datasets generated from a single man.
"""

import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Initializing parser to argparse.ArgumentParser
parser = argparse.ArgumentParser( description = __doc__)
# Specificing that the program requires the '--vcf' informer along with a VCF file
parser.add_argument('-vcf', '--vcf', help = 'Requires VCF file path')
# Initializing args to parse through all arguments
args = parser.parse_args()

# Saves the input file as a varaible
sperm_genotype_calls_file = open(args.vcf, 'r')


# Initializing lists to store all read depths for all datasets
##############################
blood_standard_dp_list = []
live_standard_dp_list = []
dead_standard_dp_list = []

live_MALBAC1_dp_list = []
live_MALBAC2_dp_list = []
live_MALBAC3_dp_list = []

dead_MALBAC1_dp_list = []
dead_MALBAC2_dp_list = []
dead_MALBAC3_dp_list = []
##############################

# Initializing variables to keep track of the number of variants with unknown coverage for each dataset
##############################
blood_standard_unknown_variants = 0
live_standard_unknown_variants = 0
dead_standard_unknown_variants = 0

live_MALBAC1_unknown_variants = 0
live_MALBAC2_unknown_variants = 0
live_MALBAC3_unknown_variants = 0

dead_MALBAC1_unknown_variants = 0
dead_MALBAC2_unknown_variants = 0
dead_MALBAC3_unknown_variants = 0
##############################

# Initializing variables to keep track of the number of allelic dropout events for each dataset
##############################
blood_standard_het_gt = 0
live_standard_dropout = 0
dead_standard_dropout = 0

live_MALBAC1_dropout = 0
live_MALBAC2_dropout = 0
live_MALBAC3_dropout = 0

dead_MALBAC1_dropout = 0
dead_MALBAC2_dropout = 0
dead_MALBAC3_dropout = 0
##############################

# Initializing variables to keep track of the number of 'live' and 'dead' variants
##############################
live_variant_count = 0
live_transition_count = 0
live_transversion_count = 0

dead_variant_count = 0
dead_transition_count = 0
dead_transversion_count = 0
##############################



# Parsing through every variant line in a given VCF to extract various amounts of information
for line in sperm_genotype_calls_file:
	# Remove newline characters at the beginning and end of the line
	line = line.rstrip()


	# Skip all lines that start with a '#' to ensure that only variant lines are read in and not the headers
	if not line.startswith("#"):
		# Split the line via whitespace to be able to call each column via index number
		line = line.split()

		# Given a line, I am isolating each datasets column respectively
		##############################
		blood_standard_column = line[9]
		live_standard_column = line[10]
		dead_standard_column = line[11]
		
		live_MALBAC1_column = line[12]
		live_MALBAC2_column = line[13]
		live_MALBAC3_column = line[14]
		
		dead_MALBAC1_column = line[15]
		dead_MALBAC2_column = line[16]
		dead_MALBAC3_column = line[17]
		##############################


		# Making a list of all read depth values and calculating the number of variants with unknown coverages for each dataset
		############################## 
		if blood_standard_column == './.' or blood_standard_column.split(':')[2] == '.':
			blood_standard_unknown_variants += 1
		else:
			blood_standard_dp_list.append(float(blood_standard_column.split(':')[2]))

		if live_standard_column == './.' or live_standard_column.split(':')[2] == '.':
			live_standard_unknown_variants += 1
		else:
			live_standard_dp_list.append(float(live_standard_column.split(':')[2]))

		if dead_standard_column == './.' or dead_standard_column.split(':')[2] == '.':
			dead_standard_unknown_variants += 1
		else:
			dead_standard_dp_list.append(float(dead_standard_column.split(':')[2]))


		if live_MALBAC1_column == './.' or live_MALBAC1_column.split(':')[2] == '.':
			live_MALBAC1_unknown_variants += 1 
		else:
			live_MALBAC1_dp_list.append(float(live_MALBAC1_column.split(':')[2]))

		if live_MALBAC2_column == './.' or live_MALBAC2_column.split(':')[2] == '.':
			live_MALBAC2_unknown_variants += 1
		else:
			live_MALBAC2_dp_list.append(float(live_MALBAC2_column.split(':')[2]))

		if live_MALBAC3_column == './.' or live_MALBAC3_column.split(':')[2] == '.':
			live_MALBAC3_unknown_variants += 1
		else:
			live_MALBAC3_dp_list.append(float(live_MALBAC3_column.split(':')[2]))


		if dead_MALBAC1_column == './.' or dead_MALBAC1_column.split(':')[2] == '.':
			dead_MALBAC1_unknown_variants += 1
		else:
			dead_MALBAC1_dp_list.append(float(dead_MALBAC1_column.split(':')[2]))

		if dead_MALBAC2_column == './.' or dead_MALBAC2_column.split(':')[2] == '.':
			dead_MALBAC2_unknown_variants += 1
		else:
			dead_MALBAC2_dp_list.append(float(dead_MALBAC2_column.split(':')[2]))

		if dead_MALBAC3_column == './.' or dead_MALBAC3_column.split(':')[2] == '.':
			dead_MALBAC3_unknown_variants += 1
		else:
			dead_MALBAC3_dp_list.append(float(dead_MALBAC3_column.split(':')[2]))
		##############################


		# Using the heterozygous sites identified in the blood standard dataset to calculate the rate of allelic dropout for each remaining dataset
		##############################
		# Initializing variable to gain access solely to the dataset's genotypes
		blood_standard_gt = blood_standard_column.split(':')[0]
		live_standard_gt = live_standard_column.split(':')[0]
		dead_standard_gt = dead_standard_column.split(':')[0]
		live_MALBAC1_gt = live_MALBAC1_column.split(':')[0]
		live_MALBAC2_gt = live_MALBAC2_column.split(':')[0]
		live_MALBAC3_gt = live_MALBAC3_column.split(':')[0]
		dead_MALBAC1_gt = dead_MALBAC1_column.split(':')[0]
		dead_MALBAC2_gt = dead_MALBAC2_column.split(':')[0]
		dead_MALBAC3_gt = dead_MALBAC3_column.split(':')[0]

		# Testing to ensure only heterozygous sites in blood standard are used for testing
		if blood_standard_gt == '0/1':
			blood_standard_het_gt += 1

			#####print(blood_standard_column.split(':')[0] + '\t' + live_standard_column.split(':')[0] + '\t' + dead_standard_column.split(':')[0] + '\t' + live_MALBAC1_column.split(':')[0] + '\t' + live_MALBAC2_column.split(':')[0] + '\t' + live_MALBAC3_column.split(':')[0] + '\t' + dead_MALBAC1_column.split(':')[0] + '\t' + dead_MALBAC2_column.split(':')[0] + '\t' + dead_MALBAC3_column.split(':')[0])

			# Testing every dataset for allelic dropout and incrementing their count if necessary
			if blood_standard_gt != live_standard_gt and live_standard_gt != './.':
				live_standard_dropout += 1

			if blood_standard_gt != dead_standard_gt and dead_standard_gt != './.':
				dead_standard_dropout += 1

			if blood_standard_gt != live_MALBAC1_gt and live_MALBAC1_gt != './.':
				live_MALBAC1_dropout += 1

			if blood_standard_gt != live_MALBAC2_gt and live_MALBAC2_gt != './.':
				live_MALBAC2_dropout += 1

			if blood_standard_gt != live_MALBAC3_gt and live_MALBAC3_gt != './.':
				live_MALBAC3_dropout += 1

			if blood_standard_gt != dead_MALBAC1_gt and dead_MALBAC1_gt != './.':
				dead_MALBAC1_dropout += 1

			if blood_standard_gt != dead_MALBAC2_gt and dead_MALBAC2_gt != './.':
				dead_MALBAC2_dropout += 1

			if blood_standard_gt != dead_MALBAC3_gt and dead_MALBAC3_gt != './.':
				dead_MALBAC3_dropout += 1
		##############################


		# Finding the set of 'live' variants and 'dead' variants
		# Calculating the number of transitions, tranversions, and the transition/transversion ratio for each set of variants
		##############################
		# Initializing list to contain the live and dead malbac replicates genotype information
		live_MALBAC_gt_list = [live_MALBAC1_gt, live_MALBAC2_gt, live_MALBAC3_gt]
		dead_MALBAC_gt_list = [dead_MALBAC1_gt, dead_MALBAC2_gt, dead_MALBAC3_gt]

		# Initializing variables to hold the reference and alternative bases reported for variants
		ref_base = line[3]
		alt_base = line[4]
		
		# Using an IF statement to determine the 'live' variants based on the parameters within the statement
		if blood_standard_gt == '0/0' and (live_MALBAC_gt_list.count('0/1') >= 2 or live_MALBAC_gt_list.count('1/1') >= 2) and  dead_MALBAC_gt_list.count('0/0') >= 2:
			live_variant_count += 1

			# Using an IF/ELSE statement to determine and count transitions, else count as transversions
			if (ref_base == 'A' and alt_base == 'G') or (ref_base == 'G' and alt_base == 'A') or (ref_base == 'C' and alt_base == 'T') or (ref_base == 'T' and alt_base == 'C'):
				live_transition_count += 1
			else:
				live_transversion_count += 1

		# Using an IF statement to determine the 'dead' variants based on the parameters within the statement
		if blood_standard_gt == '0/0' and live_MALBAC_gt_list.count('0/0') >= 2 and ( dead_MALBAC_gt_list.count('0/1') >= 2 or dead_MALBAC_gt_list.count('1/1') >= 2 ):
			dead_variant_count += 1

			# Using an IF/ELSE statement to determine and count transitions, else count as transversions
			if (ref_base == 'A' and alt_base == 'G') or (ref_base == 'G' and alt_base == 'A') or (ref_base == 'C' and alt_base == 'T') or (ref_base == 'T' and alt_base == 'C'):
				dead_transition_count += 1
			else:
				dead_transversion_count += 1


			#####print(blood_standard_column.split(':')[0] + '\t' + live_MALBAC1_column.split(':')[0] + '\t' + live_MALBAC2_column.split(':')[0] + '\t' + live_MALBAC3_column.split(':')[0] + '\t' + dead_MALBAC1_column.split(':')[0] + '\t' + dead_MALBAC2_column.split(':')[0] + '\t' + dead_MALBAC3_column.split(':')[0])
		##############################


#####print(str(blood_standard_het_gt) + '\t' + str(live_standard_dropout) + '\t' + str(dead_standard_dropout) + '\t' + str(live_MALBAC1_dropout) + '\t' + str(live_MALBAC2_dropout) + '\t' + str(live_MALBAC3_dropout) + '\t' + str(dead_MALBAC1_dropout) + '\t' + str(dead_MALBAC2_dropout) + '\t' + str(dead_MALBAC3_dropout))


# Printing to standard output the dataset's mean in variant coverage, variance in variant coverage, number of variants with unknown coverages and etc.
##############################
print('Dataset_Name\tMean\tVariance\tUnknown\tAllelic_Dropout\tAllelic_Dropout_Rate')
print('Blood Standard\t' + str(round(np.mean(blood_standard_dp_list), 2)) + '\t' + str(round(np.var(blood_standard_dp_list), 5)) + '\t' + str(blood_standard_unknown_variants) + '\t\'TRUTH\'\t\'TRUTH\'' )
print('Live Standard\t' + str(round(np.mean(live_standard_dp_list), 2)) + '\t' + str(round(np.var(live_standard_dp_list), 4)) + '\t' + str(live_standard_unknown_variants) + '\t' + str(live_standard_dropout) + '\t' + str(round(live_standard_dropout/blood_standard_het_gt, 2)))
print('Dead Standard\t' + str(round(np.mean(dead_standard_dp_list), 2)) + '\t' + str(round(np.var(dead_standard_dp_list), 3)) + '\t' + str(dead_standard_unknown_variants) + '\t' + str(dead_standard_dropout) + '\t' + str(round(dead_standard_dropout/blood_standard_het_gt, 2)))
print('Live MALBAC1\t' + str(round(np.mean(live_MALBAC1_dp_list), 2)) + '\t' + str(round(np.var(live_MALBAC1_dp_list), 4)) + '\t' + str(live_MALBAC1_unknown_variants) + '\t' + str(live_MALBAC1_dropout) + '\t' + str(round(live_MALBAC1_dropout/blood_standard_het_gt, 2)))
print('Live MALBAC2\t' + str(round(np.mean(live_MALBAC2_dp_list), 2)) + '\t' + str(round(np.var(live_MALBAC2_dp_list), 4)) + '\t' + str(live_MALBAC2_unknown_variants) + '\t' + str(live_MALBAC2_dropout) + '\t' + str(round(live_MALBAC2_dropout/blood_standard_het_gt, 2)))
print('Live MALBAC3\t' + str(round(np.mean(live_MALBAC3_dp_list), 2)) + '\t' + str(round(np.var(live_MALBAC3_dp_list), 4)) + '\t' + str(live_MALBAC3_unknown_variants) + '\t' + str(live_MALBAC3_dropout) + '\t' + str(round(live_MALBAC3_dropout/blood_standard_het_gt, 2)))
print('Dead MALBAC1\t' + str(round(np.mean(dead_MALBAC1_dp_list), 2)) + '\t' + str(round(np.var(dead_MALBAC1_dp_list), 4)) + '\t' + str(dead_MALBAC1_unknown_variants) + '\t' + str(dead_MALBAC1_dropout) + '\t' + str(round(dead_MALBAC1_dropout/blood_standard_het_gt, 2)))
print('Dead MALBAC2\t' + str(round(np.mean(dead_MALBAC2_dp_list), 2)) + '\t' + str(round(np.var(dead_MALBAC2_dp_list), 4)) + '\t' + str(dead_MALBAC2_unknown_variants) + '\t' + str(dead_MALBAC2_dropout) + '\t' + str(round(dead_MALBAC2_dropout/blood_standard_het_gt, 2)))
print('Dead MALBAC3\t' + str(round(np.mean(dead_MALBAC3_dp_list), 2)) + '\t' + str(round(np.var(dead_MALBAC3_dp_list), 4)) + '\t' + str(dead_MALBAC3_unknown_variants) + '\t' + str(dead_MALBAC3_dropout) + '\t' + str(round(dead_MALBAC3_dropout/blood_standard_het_gt, 2)))
print('\nThe number of heterozygous sites identified in the blood standard dataset (denominator): ' + str(blood_standard_het_gt))
print('\nVariant\tCount\tTransI\tTransV\tTransitions/Transversions')
print('\'Live\'\t' + str(live_variant_count) + '\t' + str(live_transition_count) + '\t' + str(live_transversion_count) + '\t' + str(round(live_transition_count/live_transversion_count, 2)))
print('\'Dead\'\t' + str(dead_variant_count) + '\t' + str(dead_transition_count) + '\t' + str(dead_transversion_count) + '\t' + str(round(dead_transition_count/dead_transversion_count, 2)))
##############################



# Creating a boxplot for each dataset to show the coverage distribution in the same figure
##############################
fig = plt.figure()# Start a new figure
# Creating a list to store all dataset read depth list
dataset = [blood_standard_dp_list, live_standard_dp_list, dead_standard_dp_list, live_MALBAC1_dp_list, live_MALBAC2_dp_list, live_MALBAC3_dp_list, dead_MALBAC1_dp_list, dead_MALBAC2_dp_list, dead_MALBAC3_dp_list]
dataset_names = ['Blood Standard', 'Live Standard', 'Dead Standard', 'Live MALBAC1', 'Live MALBAC2', 'Live MALBAC3', 'Dead MALBAC1', 'Dead MALBAC2', 'Dead MALBAC3']
plt.boxplot(dataset)# Plotting boxplots with the given data
plt.xticks(range(9), dataset_names, rotation = 17)# Labeling each dataset with it's corresponding boxplot
plt.tight_layout(pad = 2)# Creating a pad to allow the x-axis label to be viewed, Thanks to Steven :)
plt.xlabel('Datasets')# X-axis label
plt.ylabel('Read depth')# Y-axis label
plt.title('Variant Coverage Distribution')# Title of plot
fig.savefig('variant_coverage_distributions.png')# Saving figure
##############################


		







