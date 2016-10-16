#!/usr/bin/env python3
"""Takes a fasta of contigs as input, scans all six reading frames for ORFs that start with an ATG,
end with a STOP codon, and are >100bp in length (not including the stop codon), and reports the
longest ORF. Note: Since a stop codon does not code for a protein, it shouldn't count towards the
length of an ORF.

Usage: python3 call_orfs.py <contigs file>

Args:
	contigs file = Path to the contigs file.

Outputs:
	all_orfs.fna = A fasta formatted file with nucleotide sequences from predicted ORFs
	all_proteins.faa = A fasta formatted file with amino acid sequences from predicted ORFs that
have been translated into an amino acid sequence.
"""


import sys
import collections


# sys.arg is a list containing 2 elements: the script name and 1 command line argument
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 2):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)

# This function takes the reverse complement of a given sequence
def reverse_complement(dna_seq):
    complement_dictionary = {'A':'T','T':'A','G':'C','C':'G'}
    return "".join([complement_dictionary[base] for base in reversed(dna_seq)])


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



# Save the input argument as contigs_file varaible
contigs_file = open(sys.argv[1], "r")

# Initializing contigs dictionary
contigs_dictionary = {}

# Use a for loop to iterate through the file line by line
for line in contigs_file:
	# Remove newline
	line = line.rstrip()
	# (Getting the sequence name) If the line starts with the character ">" then,
	if line.startswith(">"):
		# Assign seq_name
		seq_name = line
	# (Read in sequence in uppercase) If the line does not start with the character ">" then,
	if not line.startswith(">"):
		# Make all of the characters within the line uppercase
		line = line.upper()
	# Put the seq_name and sequence in contigs_dictionary
	contigs_dictionary[seq_name] = line



# The key: seq ID and value: orf sequence
all_orfs_dictionary = {}
all_orfs_dictionary = collections.OrderedDict()
# Use a for loop to iterate through the keys and values of contigs_dictionary.items()
for contig_name,contig_sequence in contigs_dictionary.items():
	# Assigning the reverse complement of the contig_sequence to contig_reverse_complement
	contig_reverse_complement = reverse_complement(contig_sequence)
	# Creating a list of contig_sequence and contig_reverse_complement names seq
	seqs = [contig_sequence, contig_reverse_complement]
	# Assigning a sequence orientation number to keep track of when the orientation changes
	sequence_orientation = 0
	# Use a for loop to iterate through seqs list
	for seq in seqs:
		# Increment orientation (Note: 1 codes for forward and 2 codes for reverse)
		sequence_orientation += 1
		# Use a for loop to create a sliding window reading frame shifting for each frame
		for reading_frame in range(0,3):
			# Assigning currently_in_orf to False
			currently_in_orf = False
			# Creating a variable count to keep track of ORF to maintain specific naming conventions
			count = 0
			# Use a for loop to iterate through the sequence taking 3 stepsizes based on the current reading frame
			for i in range(reading_frame, len(seq), 3):
				# Use a sliding window to hold the current_codon within the sequence
				current_codon = seq[i:i+3]
				# IF the current_codon is 'ATG' then,
				if current_codon == 'ATG':
					# IF currently_in_orf is FALSE then,
					if currently_in_orf == False:
						# Set currently_in_orf to TRUE
						currently_in_orf = True
						# Assign orf_start to the current position
						orf_start = i
				# IF current_codon is a STOP codon then,
				if (current_codon == 'TAA') or (current_codon == 'TAG') or (current_codon == 'TGA'):
					# IF currently_in_orf is TRUE then,
					if currently_in_orf == True:
						# Set orf_stop to the current position i
						orf_stop = i
						# Increment count by 1
						count += 1
						# Add current ORF to all_orf_dictionary with a custom name and it's corresponding sequence
						all_orfs_dictionary[contig_name+"_orientation_"+str(sequence_orientation)+"_rf_"+str(reading_frame)+"_orf_"+str(count)] = seq[orf_start:orf_stop+3]# ADD IN UNIQUE KEY NAME TO DESIGNATE IF I AM IN THE FORWARD OR REVERSE
						# Set currently_in_orf to FALSE
						currently_in_orf = False



# Open all_orfs.fna output file for writing
all_orfs_output = open("all_orfs.fna", "w")

# Open all_proteins.faa output file for writing
all_proteins_output = open("all_proteins.faa", "w")

# Intialize an orf_count variable
orf_count = 0
# Use a for loop to iterate through all_orfs_dictionary.items()
for contig_name,contig_orf in all_orfs_dictionary.items():
	# IF ORF is greater or equal to 100bp excluding the stop codon then,
	if len(contig_orf)-3 >= 100:
		# Increment orf_count
		orf_count += 1
		# Print the current ORF to all_orf_output
		print(contig_name+'\n'+contig_orf, file = all_orfs_output)
		# Intialize an protein_sequence variable
		protein_sequence = ''
		# Use a for loop to iterate through the ORF with a stepsize of 3 exlcuding the stop codon
		for i in range(0,len(contig_orf)-3,3):
			# IF the contig_orf within the sliding window exist within the codons_to_aa dictionary then,
			if contig_orf[i:i+3] in codons_to_aa:
				# Translate the current codon to an aa and append it to the protein_sequence
				protein_sequence += codons_to_aa[contig_orf[i:i+3]]
		# Print the current protein sequence to all_proteins_output
		print(contig_name+'\n'+protein_sequence, file = all_proteins_output)
# Print the number of ORFs
print(orf_count)

# Close the output files
all_orfs_output.close()
all_proteins_output.close()
