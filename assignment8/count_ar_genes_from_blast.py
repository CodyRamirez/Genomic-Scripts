#!/usr/bin/env python3
"""Takes blast_to_card.txt as input and calculates the number of antibiotic resistance genes identified in
contigs using an identity threshold of >80% amino acid identity over >85% of the subject sequence.

Usage: python3 count_ar_genes_from_blast.py <blast output>

Args:
	blast output = Path to blast_to_card

Outputs:
	The number of antibiotic resistance genes identified in contigs using an identity threshold of >80% amino
	acid identity over >85% of the subject sequence.
"""


import sys


# sys.arg is a list containing 2 elements: the script name and 1 command line argument
# Check that all the command line argument was given. If not, prints an error statement, the documentation and then exits.
if (len(sys.argv) != 2):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)

# Save the input argument as contigs_file varaible
blast_to_card_file = open(sys.argv[1], "r")
# Intialize a count variable
count = 0
# Use a for loop to iterate through the file line by line
for line in blast_to_card_file:
	# Remove newline
	line = line.rstrip()
	# IF an identity threshold of >80% amino acid identity over >85% of the subject sequence then,
	if (float(line.split()[2]) > 80) and ((float(line.split()[3])/float(line.split()[12]))*100 > 85):
		# Print out the line
		print('\n'+line)
		# Increment count
		count += 1
# Summary print statement
print('\nThe number of antibiotic resistance genes identified in contigs using an identity threshold of >80% amino acid identity over >85% of the subject sequence is: '+str(count))