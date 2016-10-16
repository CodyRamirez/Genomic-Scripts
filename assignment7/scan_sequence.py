#!/usr/bin/env python3
"""Scan a DNA sequence to find putative binding sites

Usage: python3 scan_sequence.py <scoring_matrix> <sequence_file> <score_threshold>

Args:
	scoring_matrix = Path to scoring matrix. The rows of the matrix correspond to 
		A, C, G, and T, and the columns correspond to positions
	sequence_file = Path to DNA sequence file.
	score threshold = An integer
"""
import sys

# Helper dictionary for reverse complementation
reverse_comp = { 'A':'T', 'C':'G', 'G':'C', 'T':'A' }

###############################################################
# Begin functions
###############################################################

def create_scoring_matrix_from_file(matrix_file):
	"""Creates a scoring matrix given a matrix file. The scoring matrix is a list of dictionaries. 
		The index of the list corresponds to the column of the matrix file.
		Each key,value pair in the dictionary is the base which corresponds to a row from the matrix file
		The values assigned to each key is the value within the matrix file corresponding to that row (base).

	Args:
		matrix file: The columns are positions and the rows correspond to A, C, G, and T.

	Returns:
		scoring_matrix: List of dictionaries where the index of the list corresponds to positions
		and the dictionaries key,value pair are base,value corresponding to the matrix rows within the column
	"""
	file_data = [ x.split() for x in open(matrix_file).readlines() ]
	scoring_matrix = [ dict(A=float(a_score), C=float(c_score), G=float(g_score), 
		T=float(t_score)) for a_score, c_score, g_score, t_score in 
		zip(file_data[0], file_data[1], file_data[2], file_data[3])
	]

	return scoring_matrix


def read_sequence_from_file(file_name):
	"""Opens a given file, reads in the first line, removes all whitespace and then returns that line.
		The first line corresponds to a given sequence.

	Args:
		file_name: A file containing a sequence

	Returns: 
		A sequence
	"""
	return open(file_name).readlines()[0].strip()



def score_with_matrix(subseq, matrix):
	"""Calculates the score of the current sequence within the whole sequence based on the given scoring matrix.

	Args:
		subseq: a subsection of the whole sequence based on the number of positions within the given matrix
		matrix: a scoring matrix; columns are positions and the rows corresponds to A, C, G, and T
	Returns:
		The score of the current sequence motif
	"""
	return sum([ score[ base ] for score, base in zip(matrix, subseq)])


def get_reverse_complement(in_string):
	"""Returns the reverse complement of a sequence

	Args:
		in_string: a string of sequence

	Returns:
		reverse complement
	"""
	return (''.join([ reverse_comp[x] for x in in_string[::-1] ]))

###############################################################
# End functions
###############################################################

###############################################################
# Begin main script
###############################################################

# Check the correct number of command line arguments
if(len(sys.argv)!= 4):
	sys.exit(__doc__)

# Assigning arguments to variables
score_matrix_file = sys.argv[1]
sequence_file = sys.argv[2]
score_threshold = float(sys.argv[3])

score_matrix = create_scoring_matrix_from_file(score_matrix_file)# Converts the score_matrix_file into a useable score_matrix
motif_width = len(score_matrix)# Assigning a motif_width based on the length of the score_matrix
search_sequence = read_sequence_from_file(sequence_file)# Assigning search_sequence to the current sequence (forward position) within a sequence_file
search_reverse_complement = get_reverse_complement(search_sequence)# Assigning seatch_reverse_complement to the current sequence (reverse position) within a sequence_file

# Calculate the number of matrix 'windows' for calculating sequence scores
last_index = len(search_sequence) - motif_width + 1

# Use a list comprehension to iterate through the sequence using a sliding window and create a triplet entry with:
# Index[0] = the sequence of the putative binding site
# Index[1] = the leftmost position of the putative binding site
# Index[2] = the score of the putative binding site
forward_hit_list_raw = [ (search_sequence[i:i+motif_width], i, score_with_matrix(search_sequence[i:i+motif_width], score_matrix)) for i in range(last_index)]

# Use a list comprehension to iterate through the foward_hit_list_raw and filter based on the threshold
forward_hit_list = [ forward_hit_list_raw[i] for i in range(len(forward_hit_list_raw)) if forward_hit_list_raw[i][2] >= score_threshold]

# Use a list comprehension to iterate through the reverse complement sequence using a sliding window and create a triplet entry with:
# Index[0] = the sequence of the putative binding site
# Index[1] = the leftmost position of the putative binding site
# Index[2] = the score of the putative binding site
reverse_hit_list_raw = [ (search_sequence[i:i+motif_width], len(search_sequence) - i - 1, score_with_matrix(search_reverse_complement[i:i+motif_width], score_matrix)) for i in range(last_index)]

# Use a list comprehension to iterate through the reverse_hit_list_raw and filter based on the threshold
reverse_hit_list = [ reverse_hit_list_raw[i] for i in range(len(reverse_hit_list_raw)) if reverse_hit_list_raw[i][2] >= score_threshold]

# IF forward_hit_list is empty then,
if len(forward_hit_list) == 0:
	print("No threshold-exceeding hits found in the forward direction!")# Print no hits found
else:# ELSE
	print("orientation\tsequecne\tposition\tscore")# Print a header containing: orientation, sequence, position, and score
	for hit in forward_hit_list:# Use a for loop to iterate through forward_hit_list
		print("forward\t\t{bind_site}\t{position:d}\t\t{score:.2f}".format(bind_site=hit[0], position=hit[1], score=hit[2]))# Print out the orientation, sequence, position and score rounded to two decimal places

# IF forward_hit_list is empty then,
if len(reverse_hit_list) == 0:
	print("No threshold-exceeding hits found in the reverse direction!")# Print not hits found
else:# ELSE
	print("orientation\tsequence\tposition\tscore")# Printer a header containing: orientation, sequence, position, and score
	for hit in reverse_hit_list:# Use a for loop to iterate through reverse_hit_list
		print("reverse\t\t{bind_site}\t{position:d}\t\t{score:.2f}".format(bind_site=hit[0], position=hit[1], score=hit[2]))# Print out the orientation, sequence, postion and score rounded to two decimal places
