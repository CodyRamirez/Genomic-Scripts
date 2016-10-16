#!/usr/bin/env python3
"""Filter variants such that you only keep records in which all three individuals
	in the trio have GQ scores at or above a given threshold. If any of the three
	are missing GQ values or if any value is below the threshold, that record 
	should not be kept. If the user does not provide a threshold, the script should
	not filter the variants.

	Before Modification: Counts the number of variants that clearly violate the 
	rules of Mendelian segregation, given the trio's relationship to one another.

	Usage: 
		python3 violate_MS.py <SNV_indel VCF> [GQ threshold*]
		*default: no thresholding based on genotype quality scores

	Args:
		SNV_indel VCF = VCF file of small genome variation
		GQ threshold = An integer threshold to filter GQ

	Outputs:
		Variants such that you only keep records in which all three individuals
	in the trio have GQ scores at or above a given threshold.

		Before Modification: The number of variants that clearly violate the 
		rules of Mendelian segregation
"""
import sys


# Intialize gq_threshold to nothing
gq_threshold = " "

# Use a for loop to iterate through the number of arguments inputed
for arg in sys.argv[1:]:
	if arg.isdigit() == True:# IF arg is a digit then,
		gq_threshold = int(arg)# Assign arg to variable gq_threshold
	else:# ELSE
		snv_indel_vcf_file = open(arg, "r")# Assign arg to snv_indel_vcf_file


################### MAIN ####################
# Intializing count variables
violate_ms_count = 0

if gq_threshold == " ":
	# Use a for loop to iterate through the file line by line to read a single line at a time
	for line in snv_indel_vcf_file:
		# Remove newline characters at the beginning and end of the line
		line = line.rstrip()
		# Skip all lines that start with a '#' to ensure that only variant lines are read in and not the headers
		if not line.startswith("#"):
			# Split the line via whitespace to be able to call each column via index number
			line = line.split()

			# Filtering out sex chromosomes to only count violates occuring within the autosomes 
			if line[0] != 'X':
				if line[0] != 'Y':

					# Access the genotype info within NA12878's column
					child_genotype_info = line[10].split(':')[0].split('/')
					
					# Accessing the genotypes info of the parents of NA12878 and merging them into a single list
					parents_genotype_info = line[23].split(':')[0].split('/') + line[24].split(':')[0].split('/')

					# Testing to see if every allele in the child exist within the pool of alleles of the parents, then
					if all([g in parents_genotype_info for g in child_genotype_info]) == True:
						# Testing cases where the number of alleles represented within the child are more than the parents pooled alleles
						if child_genotype_info.count('1') > parents_genotype_info.count('1'):
							# If so, increment violate_ms_count
							violate_ms_count += 1
						# Testing cases where the number of alleles represented within the child are more than the parents pooled alleles
						if child_genotype_info.count('0') > parents_genotype_info.count('0'):
							# If so, increment violate_ms_count
							violate_ms_count += 1
					# Testing to see if every allele in the child exist within the pool of alleles of the parents, then
					if all([g in parents_genotype_info for g in child_genotype_info]) == False:
						# If not so, then increment violate_ms_count
						violate_ms_count += 1

	print('Number of variants that clearly violate the rules of Mendelian segregation: ' + str(violate_ms_count))



else:
	# Use a for loop to iterate through the file line by line to read a single line at a time
	for line in snv_indel_vcf_file:
		# Remove newline characters at the beginning and end of the line
		line = line.rstrip()
		# Skip all lines that start with a '#' to ensure that only variant lines are read in and not the headers
		if not line.startswith("#"):
			# Split the line via whitespace to be able to call each column via index number
			line = line.split()

			# Filtering out sex chromosomes to only count violates occuring within the autosomes 
			if line[0] != 'X':
				if line[0] != 'Y':

					# Access the genotype info within NA12878's column
					child_genotype_info = line[10].split(':')[0].split('/')
					
					# Accessing the genotypes info of the parents of NA12878 and merging them into a single list
					parents_genotype_info = line[23].split(':')[0].split('/') + line[24].split(':')[0].split('/')
					# Intializing a gq_index to keep track of the GQ position
					gq_index = 0

					# Split the format column by ':' to gain access to the individual access to each variable
					format_column = line[8].split(':')
					na12878_column = line[10].split(':')
					na12891_column = line[23].split(':')
					na12892_column = line[24].split(':')

					# Ensuring that the GQ column exist
					if len(na12878_column) >= 4:
						if len(na12891_column) >= 4:
							if len(na12892_column) >= 4:

								# Using a for loop to iterate through the format column
								for i in range(0, len(format_column)):
									# If the current format column is GQ then,
									if format_column[i] == 'GQ':
										# Store the index number as gt_index to later use in the patient's field to find their GT information
										gq_index = i

										#print(format_column, na12878_column, na12891_column, na12892_column)
										# Filtering out the missing GQ values
										if na12878_column[gq_index] != '.':
											if na12891_column[gq_index] != '.':
												if na12892_column[gq_index] != '.':

													# Ensuring that each individual is >= the threshold set
													if int(na12878_column[gq_index]) >= gq_threshold:
														if int(na12891_column[gq_index]) >= gq_threshold:
															if int(na12892_column[gq_index]) >= gq_threshold:											

																# Testing to see if every allele in the child exist within the pool of alleles of the parents, then
																if all([g in parents_genotype_info for g in child_genotype_info]) == True:
																	# Testing cases where the number of alleles represented within the child are more than the parents pooled alleles
																	if child_genotype_info.count('1') > parents_genotype_info.count('1'):
																		# If so, increment violate_ms_count
																		violate_ms_count += 1
																	# Testing cases where the number of alleles represented within the child are more than the parents pooled alleles
																	if child_genotype_info.count('0') > parents_genotype_info.count('0'):
																		# If so, increment violate_ms_count
																		violate_ms_count += 1
																# Testing to see if every allele in the child exist within the pool of alleles of the parents, then
																if all([g in parents_genotype_info for g in child_genotype_info]) == False:
																	# If not so, then increment violate_ms_count
																	violate_ms_count += 1
	print('Number of variants that clearly violate the rules of Mendelian segregation: ' + str(violate_ms_count))

