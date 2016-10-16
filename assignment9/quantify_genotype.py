#!/usr/bin/env python3
"""Quantifies the number of homozygous and heterozygous SNVs and indels, for the
	individual NA12878. Count the total number of genotype calls that are:
	1) homozygous reference
	2) heterozygous
	3) homozygous alternate
	4) one or more alleles missing


	Usage: 
		python3 quantify_genotype.py <SNV_indel VCF>

	Args:
		SNV_indel VCF = VCF file of small genome variation

	Outputs:
		A table with the four possibilities and their corresponding counts.
"""
import sys

# sys.arg is a list containing 2 elements: the script name and a command line argument
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 2):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)



################### MAIN ####################
# Intializing count variables
homo_ref_genotype_count = 0
hetero_genotype_count = 0
homo_alt_genotype_count = 0
missing_genotype_count = 0

# Save the input arguments as varaibles
snv_indel_vcf_file = open(sys.argv[1], "r")

# Use a for loop to iterate through the file line by line to read a single line at a time
for line in snv_indel_vcf_file:
	# Remove newline characters at the beginning and end of the line
	line = line.rstrip()
	# Skip all lines that start with a '#' to ensure that only variant lines are read in and not the headers
	if not line.startswith("#"):
		# Split the line via whitespace to be able to call each column via index number
		line = line.split()
		# Access the genotype info within NA12878's column
		genotype_info = line[10].split(':')[0]
		
		# IF one or more alleles are missing then,
		if '.' in genotype_info:
			# Increment the missing genotype count
			missing_genotype_count += 1
		# IF genotype is homozygous reference then,
		if genotype_info == '0/0':
			# Increment the homo ref genotype
			homo_ref_genotype_count += 1
		# IF genotype is homozygoes alternative then,
		if genotype_info == '1/1':
			# Increment the homo alt genotype
			homo_alt_genotype_count += 1
		# IF genotype is heterozygous then,
		if genotype_info in ('0/1', '1/0'):
			# Increment hetero genotype
			hetero_genotype_count += 1


print('Homozygous reference count:\t' + str(homo_ref_genotype_count))
print('Heterozygous count:\t' + str(hetero_genotype_count))
print('Homozygous alternate count: ' + str(homo_alt_genotype_count))
print('Missing alleles count: ' + str(missing_genotype_count))



