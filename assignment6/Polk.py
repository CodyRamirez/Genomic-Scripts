#!/usr/bin/env python3
"""
Script that takes in amino acid sequence as input and outputs all possible DNA sequences that could encode this AA.
Usage: python3 Polk.py <peptide sequence(s)> <melting temperature>
"""

# Importing necessary modules
import sys

# Initializing a list for the number of arguments inputed
arg_list = []

# Use a for loop to iterate through the number of arguments inputed
for arg in sys.argv[1:]:
	if arg.isalpha() == True:# IF arg is only letters then
		arg_list.append(arg)# Append arg to arg_list
	else:# ELSE
		temp = float(arg)# Assign arg to the variable temp


#Standard table of codons - a dictionary of single-letter
#amino acid code to a list of codon choices
aa_to_codons = {}
aa_to_codons[ "A" ] = ["GCA", "GCC", "GCG", "GCT" ]
aa_to_codons[ "C" ] = ["TGC", "TGT" ]
aa_to_codons[ "D" ] = ["GAC", "GAT" ]
aa_to_codons[ "E" ] = ["GAA", "GAG" ]
aa_to_codons[ "F" ] = ["TTC", "TTT" ]
aa_to_codons[ "G" ] = ["GGA", "GGC", "GGG", "GGT" ]
aa_to_codons[ "H" ] = ["CAC", "CAT" ]
aa_to_codons[ "I" ] = ["ATA", "ATC", "ATT" ]
aa_to_codons[ "K" ] = ["AAA", "AAG" ]
aa_to_codons[ "L" ] = ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG" ]
aa_to_codons[ "M" ] = ["ATG" ]
aa_to_codons[ "N" ] = ["AAC", "AAT" ]
aa_to_codons[ "P" ] = ["CCA", "CCC", "CCG", "CCT" ]
aa_to_codons[ "Q" ] = ["CAA", "CAG" ]
aa_to_codons[ "R" ] = ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT" ]
aa_to_codons[ "S" ] = ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT" ]
aa_to_codons[ "T" ] = ["ACA", "ACC", "ACG", "ACT" ]
aa_to_codons[ "V" ] = ["GTA", "GTC", "GTG", "GTT" ]
aa_to_codons[ "W" ] = ["TGG" ]
aa_to_codons[ "Y" ] = ["TAC", "TAT" ]
aa_to_codons[ "*" ] = ["TAA", "TAG", "TGA" ]


# Inputs: Sequence of dna and amino acids
# Outputs: All possible variations of dna sequence that could encode the corresponding amino acid sequence given
def check_combinations(dna_string, aa_string, temp):
	# if this code is confusing to you, uncomment the print statement
	#print("Input DNA is:",dna_string,"Remaining AAs are:", aa_string, sep='\t')

	if (len(aa_string) == 0):# If the current length of the amino acid sequence is zero then,
		#Conduct all filtering steps within this block of code and then print all sequences that pass
		g_count = dna_string.count('G')# Initializing the variable g_count to track the number of Gs in the dna_string
		c_count = dna_string.count('C')# Initializing the variable c_count to track the number of Cs in the dna_string
		dna_temp = (64.9 + (41.0* (g_count + c_count - 16.4)/ len(dna_string) ) )# Initializing the variable dna_temp to calculate the melting temperature of dna_string

		#IF statement for filtering sequences out of specified tempature range
		if dna_temp >= temp-0.5 and dna_temp <= temp+0.5:# If dna_temp is within 0.5 degrees of specified temp then continue

			#IF statement for filtering sequence with restriction sites NdeI, XhoI, TaqI and BfaI
			if ('CATATG' and 'CTCGAG' and 'TCGA' and 'CTAG') not in dna_string:# If none of these sequences exist within dna_string then continue

				print(dna_string + '\t' + str(dna_temp) )# Print the dna sequence
	
	# Since aa_string still contains some aa continue translating them into dna sequence
	# Otherwise, 
	else:

		# Assigning the first index of aa_string to current_AA
		current_AA = aa_string[0];

		# Use a for loop to iterate through all values corresponding to the current_AA in aa_to_codons dictionary 
		for single_codon in aa_to_codons[current_AA]:

			# Assign new_dna_string to dna_string with the addition of single_codon
			new_dna_string = dna_string + single_codon

			# Calling on the function itself (recursively) but using the next index within the aa_string to calculate the next aa
			check_combinations(new_dna_string, aa_string[1:], temp)

### Main Script ###

for oligo in arg_list:# Use a for loop to iterate through arg_list
	print('The amino acid sequence and melting temperature used.\tOligo: ' + oligo + '\tMT: ' + str(temp))# Print out the current oligo and it's corresponding temperature
	print('DNA sequence\t\tMelting Temperature')# Print out the header for the sequences and temperatures
	check_combinations( "", oligo, temp )# Call on check_combinations to convert the current aa sequence to dna sequence and it's corresponding melting temperature
