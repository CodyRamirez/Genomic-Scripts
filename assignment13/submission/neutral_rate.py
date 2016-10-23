#!/usr/bin/env python3
""" Calculates the fraction of wobble positions that are conserved across all 4 species starting from
a conserved start codon until the end of the sequence (which is a stop codon in all 4 species). 
Consider only third position wobbles.

	Usage: 
		python3 neutral_rate.py <clustalw alignment file> <min bp conserved> <window length>

		Example:
			python3 neutral_rate.py PRE1.aln 

	Output file name:
		S_cer_conserved.txt
"""


import sys
import numpy as np
from scipy.stats import binom


# sys.arg is a list containing 2 elements: the script name and command line argument
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 4):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)


pre1_alignment_file = open(sys.argv[1], 'r')
min_bp_conserved = int(sys.argv[2])
window_length = int(sys.argv[3])


# Standard table of amino acids - a dictionary
# Codon choices code to a list amino acids
codons_to_aa = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
				'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
				'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
				'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
				'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
				'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
				'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
				'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
				'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
				'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
				'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
				'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
				'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
				'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
				'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
				'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}


skud_seq = ''
smik_seq = ''
scer_seq = ''
sbay_seq = ''
max_seq_length = 0

# Using a for loop to iterate through every line in pre1_alignment_file
for line in pre1_alignment_file:

	# Making each sequence whole after alignment
	##############################
	if line.startswith('Skud'):
		skud_seq += line.strip().split()[1]

	if line.startswith('Smik'):
		smik_seq += line.strip().split()[1]

	if line.startswith('Scer'):
		scer_seq += line.strip().split()[1]

	if line.startswith('Sbay'):
		sbay_seq += line.strip().split()[1]
	##############################

# Searching through every reading frame in a 3 bp sliding window to identify a conserved start codon until the end of the sequence (which is a stop codon in all 4 species) taking into account ORF parameters
##############################
# Iterating through reading frames
for reading_frame in range(0,3):

	#print('Reading Frame: ' + str(reading_frame) )

	currently_in_orf = False

	# Iterating through the sequence given the current reading frame
	for i in range(reading_frame,len(skud_seq),3):

		# Ensuring that the 'codon' sliding window sequence actually is a codon before checking for the conserved start ('M') and stop ('*') codons
		if (skud_seq[i:i+3] in codons_to_aa) and (smik_seq[i:i+3] in codons_to_aa) and (scer_seq[i:i+3] in codons_to_aa) and (sbay_seq[i:i+3] in codons_to_aa):

			# Converting all iviable codon sequences into its corresponding amino acid
			skud_aa = codons_to_aa[skud_seq[i:i+3]]
			smik_aa = codons_to_aa[smik_seq[i:i+3]]
			scer_aa = codons_to_aa[scer_seq[i:i+3]]
			sbay_aa = codons_to_aa[sbay_seq[i:i+3]]

			# Checking that all sequences share a conserved start codon
			if 'M' == skud_aa == smik_aa == scer_aa == sbay_aa:

				#print(str(i) + '\t' + skud_aa)

				if currently_in_orf == False:

					currently_in_orf = True

					start_position = i
					print(i)
			# Checking that all sequences share a conserved stop codon
			if '*' == skud_aa == smik_aa == scer_aa == sbay_aa:

				#print(str(i) + '\t' + skud_aa)

				if currently_in_orf == True:

					stop_position = (i + 3)

					# Making sure to only save the longest ORF (which also happens to be the ORF with the earliest conserved start codon and last conserved stop codon in the alignments)
					if (stop_position - start_position) > max_seq_length:

						#print(skud_seq[start_position:stop_position])

						max_seq_length = (stop_position - start_position)

						orf_start = start_position

						orf_stop = stop_position

					currently_in_orf = False
##############################

# Identifying the fraction of wobble positions that are conserved within the largest ORF
##############################
possible_wobble_sites = 0
conserved_sites = 0

for i in range(orf_start, orf_stop, 3):

	skud_codon = skud_seq[i:i+3]
	smik_codon = smik_seq[i:i+3]
	scer_codon = scer_seq[i:i+3]
	sbay_codon = sbay_seq[i:i+3]

	# Checking to make sure that all of the codons code for the same amino acid, but are not amino acid with a single codon key
	if codons_to_aa[skud_codon] == codons_to_aa[smik_codon] == codons_to_aa[scer_codon] == codons_to_aa[sbay_codon] and codons_to_aa[skud_codon] != 'M' and codons_to_aa[skud_codon] != 'W':

		#print('\n' + codons_to_aa[skud_codon] + '\t' + codons_to_aa[smik_codon] + '\t' + codons_to_aa[scer_codon] + '\t' + codons_to_aa[sbay_codon])

		# Checking to ensure the codons are considered a possible wobble position defined as: They code for the same amino acid AND share the first two bases in common
		if skud_codon[0:2] == smik_codon[0:2] == scer_codon[0:2] == sbay_codon[0:2]:

			#print(skud_codon[0:2] + '\t' + smik_codon[0:2] + '\t' +  scer_codon[0:2] + '\t' + sbay_codon[0:2])

			possible_wobble_sites += 1

			# Checking to ensure that the codons are conserved if they are all exactly the same sequence
			if skud_codon == smik_codon == scer_codon == sbay_codon:

				conserved_sites += 1
##############################

#print(max_seq_length)
#print(max_seq_length/3)
#print(skud_seq[orf_start:orf_stop])

print('Number of conserved sites: ' + '\t' + str(conserved_sites))
print('Number of wobble sites: ' + '\t' + str(possible_wobble_sites))
print('Fraction of conserved wobble sites: ' + '\t' + str(round(conserved_sites/possible_wobble_sites, 2)))



#for i in range(0,11):
#	print( str(i) + '\t' + str(binom.pmf(i, 10, 0.45)))



promoter_position_list = []

# Using an overlapping windows approach to identify promoter sequence regions >= 10 bp within each 10 sub-sequence is more conserved than expected
##############################
# Using an overlapping windows approach to identify a list of the beginning promoter prositions 
for i in range(0, orf_start - window_length):

	base_conservation_count = 0

	for base in range(i, i + window_length):
		
		if skud_seq[base] == smik_seq[base] == scer_seq[base] == sbay_seq[base]:

			base_conservation_count += 1

	if base_conservation_count >= min_bp_conserved:

		promoter_position_list.append(i)

		#print(str(i) + '\t' + skud_seq[i:i+10])



# Using the promoter position list to determine promoter sequence regions >= 10 bp within each 10 sub-sequence is more conserved than expected
output_fileobject = open('S_cer_conserved.txt', 'w')
print('Start\tStop\tScer Sequence', file = output_fileobject)

currently_in_promoter = False
complete_promoter = False

for i in range(len(promoter_position_list)-1):

	if promoter_position_list[i+1] - promoter_position_list[i] == 1:
	
		if currently_in_promoter == False:

			currently_in_promoter = True

			complete_promoter = False

			promoter_start = promoter_position_list[i]

	else:

		if currently_in_promoter == True:

			promter_stop = promoter_position_list[i]+10

			currently_in_promoter = False

			complete_promoter = True

		print(str(promoter_start) + '\t' + str(promter_stop) + '\t' + scer_seq[promoter_start:promter_stop], file = output_fileobject)

	if complete_promoter == False and i == ( len(promoter_position_list) - 2):

		promter_stop = promoter_position_list[-1] + 10

		print(str(promoter_start) + '\t' + str(promter_stop) + '\t' + scer_seq[promoter_start:promter_stop], file = output_fileobject)

output_fileobject.close()
##############################








