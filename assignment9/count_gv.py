#!/usr/bin/env python3
"""Counts the different classes of genome variation (SNVs, small indels, as well
	as larger structural variants) for the individual NA12878. The script will 
	read in both files and output a table with the names and numbers of non-reference
	alleles belonging to different variant classes, i.e. SNVs, indels, large
	deletions, tandem duplications, inversions, MEIs and BNDs.

	Usage: 
		python3 count_gv.py <SNV_indel VCF> <SV VCF>

	Args:
		SNV_indel VCF = VCF file of small genome variation
		SV VCF = VCF file of structural variation

	Outputs:
		A table with the names and numbers of non-reference alleles belonging to
		different variant classes, i.e. SNVs, indels, large deletions, tandem 
		duplications, inversions, MEIs and BNDs.
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
snv_indel_vcf_file = open(sys.argv[1], "r")
sv_vcf_file = open(sys.argv[2], "r")


#################### FUNCTIONS ####################
def non_reference_allele_counter(line):
	"""Counts the number of non-reference alleles and grabs the length of the of indels

	Args:
		line: The current line within a VCF file

	Returns:
		count
		length
	"""
	# Intializing gt_index to maintain track of the GT index within the format column(8) to find in the patient information
	count = 0
	length = 0
	gt_index = 0
	# Split the format column by ':' to gain access to the individual access to each variable
	format_column = line[8].split(':')
	# Using a for loop to iterate through the format column
	for i in range(0, len(format_column)):
		# If the current format column is GT then,
		if format_column[i] == 'GT':
			# Store the index number as gt_index to later use in the patient's field to find their GT information
			gt_index = i
			# Split na12878 column by ':' to gain access to individual index groups
			na_12878_column = line[10].split(':')
			# IF a non-reference allele exist then,
			if '1' in na_12878_column[gt_index]:
				# Calculate the length of the indel
				length = abs(len(line[3]) - len(line[4]))
				# Count the number of non-reference alleles
				count = na_12878_column[gt_index].count('1')
				
	return (count,length)



def grab_sv_length(line):
	"""Grabs the length of the sv

	Args:
		line: The current line within a VCF file

	Returns:
		length
	"""
	# Intializing gt_index to maintain track of the GT index within the format column(8) to find in the patient information
	length = 0
	gt_index = 0
	# Split the format column by ':' to gain access to the individual access to each variable
	format_column = line[8].split(':')
	# Using a for loop to iterate through the format column
	for i in range(0, len(format_column)):
		# If the current format column is GT then,
		if format_column[i] == 'GT':
			# Store the index number as gt_index to later use in the patient's field to find their GT information
			gt_index = i
			# Split na12878 column by ':' to gain access to individual index groups
			na_12878_column = line[10].split(':')
			# IF a non-reference allele exist then,
			if '1' in na_12878_column[gt_index]:
				# IF the SVLEN is a digit then,
				if line[7].split(';')[2].split('=')[1].isdigit() == True:
					# Split the INFO column by ';', grab the third element (SVLEN column), then split by '=' to access the SVLEN value and assign it to length
					length = int(line[7].split(';')[2].split('=')[1])
				# IF the SVLEN is NOT a digit then,
				if line[7].split(';')[2].split('=')[1].isdigit() == False:
					# Split the INFO column by ';', grab the third element (SVLEN column), then split by '=' to access the SVLEN value and assign its absolute value to length
					length = abs(int(line[7].split(';')[2].split('=')[1]))

	return length



#################### MAIN ####################
# Intializing count variables
snv_count = 0
indel_count = 0
del_count = 0
dup_count = 0
inv_count = 0
mei_count = 0
bnd_count = 0

# Intializing total variables
snv_total = 0
indel_total = 0
del_total = 0
dup_total = 0
inv_total = 0
mei_total = 0
bnd_total = 0

# Intializing a list to contain lengths of SV
indel_lengths_list = []
del_lengths_list = []
mei_lengths_list = []


# Use a for loop to iterate through the file line by line to read a single line at a time
for line in snv_indel_vcf_file:
	# Remove newline characters at the beginning and end of the line
	line = line.rstrip()
	# Skip all lines that start with a '#' to ensure that only variant lines are read in and not the headers
	if not line.startswith("#"):
		# Split the line via whitespace to be able to call each column via index number
		line = line.split()
		
		# IF the REF length is equal to the ALT length (indication of a SNV) then,
		if len(line[3]) == len(line[4]):
			# Count all non-reference alleles
			snv_count, length = non_reference_allele_counter(line)
			snv_total += snv_count

		# IF the REF length does not equal the ALT length (indication of an INDEL) then,
		if (len(line[3]) != len(line[4])):
			# Count all non-reference alleles
			indel_count, length = non_reference_allele_counter(line)
			indel_total += indel_count
			indel_lengths_list.append(length)
						


# Use a for loop to iterate through the file line by line to read a single line at a time
for line in sv_vcf_file:
	# Remove newline characters at the beginning and end of the line
	line = line.rstrip()
	# Skip all lines that start with a '#' to ensure that only variant lines are read in and not the headers
	if not line.startswith("#"):
		# Split the line via whitespace to be able to call each column via index number
		line = line.split()

		# IF the first element within the info column is SVTYPE=DEL then,
		if line[7].split(';')[0] == 'SVTYPE=DEL':
			# Count all non-reference alleles
			del_count, length = non_reference_allele_counter(line)
			# Incrementing the total count
			del_total += del_count

			# Assigning del_length to the current length of the DEL
			del_length = grab_sv_length(line)
			# Appending the current length to the list of DEL lengths
			del_lengths_list.append(del_length)

		# IF the first element within the info column is SVTYPE=DUP then,
		if line[7].split(';')[0] == 'SVTYPE=DUP':
			# Count all non-reference alleles
			dup_count, length = non_reference_allele_counter(line)
			# Incrementing the total count
			dup_total += dup_count

		# IF the first element within the info column is SVTYPE=INV then,
		if line[7].split(';')[0] == 'SVTYPE=INV':
			# Count all non-reference alleles
			inv_count, length = non_reference_allele_counter(line)
			# Incrementing the total count
			inv_total += inv_count

		# IF the first element within the info column is SVTYPE=MEI then,
		if line[7].split(';')[0] == 'SVTYPE=MEI':
			# Count all non-reference alleles
			mei_count, length = non_reference_allele_counter(line)
			# Incrementing the total count
			mei_total += mei_count

			# Assigning mei_length to the current length of the MEI
			mei_length = grab_sv_length(line)
			# Appending the current length to the list of MEI lengths
			mei_lengths_list.append(mei_length)

		# IF the first element within the info column is SVTYPE=BND then,
		if line[7].split(';')[0] == 'SVTYPE=BND':
			# Count all non-reference alleles
			bnd_count, length = non_reference_allele_counter(line)
			# Incrementing the total count
			bnd_total += bnd_count


# Calculating the total number of variants reported
total_count = snv_total + indel_total + del_total + dup_total + inv_total + mei_total + bnd_total
# Calculating the proportion of genomic variants that are SVs
proportion_sv = ( (del_total + dup_total + inv_total + mei_total + bnd_total)/total_count)


# Printing out a table with the names and numbers of non-reference alleles belonging to different variants classes
print('\nVariants:\tCount')
print('SNVs\t\t' + str(snv_total))
print('INDELs\t\t' + str(indel_total))
print('DELs\t\t' + str(del_total))
print('DUPs\t\t' + str(dup_total))
print('INVs\t\t' + str(inv_total))
print('MEIs\t\t' + str(mei_total))
print('BNDs\t\t' + str(bnd_total/2))
# Printing out the proportion of genomic variants that are SVs
print('\nProportion of genomic variants that are SVs: ' + str(proportion_sv))


# Creating a histogram of small INDELs, size defined/restricted from 1 to 50
fig = plt.figure()# Start a new figure
plt.hist(indel_lengths_list, bins = 100, range = (1,50), histtype = 'bar', align = 'mid')# Plotting the length distribution of small INDELs using 100 bins with bars aligned to the middle
plt.xlabel('Lengths')# X-axis label
plt.ylabel('Counts')# Y-axis label
plt.title('Length Distribution of Small INDELs')# Title of plot
fig.savefig('histogram_indels.png')# Saving figure

# Creating a histogram of large deletions DELs, size defined/restricted from 50 to 15,000
fig = plt.figure()# Start a new figure
plt.hist(del_lengths_list, bins = 1000, range = (50, 15000), histtype = 'bar', align = 'mid')# Plotting the length distribution of DELs using 1000 bins with bars aligned to the middle
plt.xlabel('Lengths')# X-axis label
plt.ylabel('Counts')# Y-axis label
plt.title('Length Distribution of DELs')# Title of plot
fig.savefig('histogram_deletions.png')# Saving figure

# Creating a histogram of mobile element insertions MEIs, size defined/restricted from 1 to 7,000
fig = plt.figure()# Start a new figure
plt.hist(mei_lengths_list, bins = 1000, range = (1, 7000), histtype = 'bar', align = 'mid')# Plotting the length distribution of MEIs using 1000 bins with bars aligned to the middle
plt.xlabel('Lengths')# X-axis label
plt.ylabel('Counts')# Y-axis label
plt.title('Length Distribution of MEIs')# Title of plot
fig.savefig('histogram_meis.png')# Saving figure
