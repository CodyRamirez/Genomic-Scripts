#!/usr/bin/env python3

""" 
nuc_count.py counts nucleotides in a fasta file

Usage: python3 nuc_count.py <fasta>

<fasta> = path to a fasta file
""" 

# Import modules
import sys
from decimal import Decimal
import collections
import re

# sys.arg is a list containing 2 elements: the script name and 1 command line argument
# Check that all the command line argument was given. If not, print the documentation and exit.
if (len(sys.argv) != 2):
    sys.exit(__doc__) 

# Save the input arguments as variables
fasta = sys.argv[1]

# Initialize a nucleotide string
nucleotides = ""

# Load the fasta sequence
# NOTE: this script assumes there is only *one* sequence in the fasta file
# Open the fasta file
with open(fasta) as f:
    # For each line in the file
    for line in f:
        # Skip lines starting with ">"
        if not line.startswith(">"):
            # Add each line to the nucleotide string
            nucleotides += line.rstrip()

# Make the nucleotide string all capital letters 
nucleotides = nucleotides.upper()

# Count the nucleotides and print output
num_a = nucleotides.count('A')
num_c = nucleotides.count('C')
num_g = nucleotides.count('G')
num_t = nucleotides.count('T')
num_n = nucleotides.count('N')

#Total of all nucleotides
total_nuc = num_a + num_c + num_g + num_t

"""
#Seeing all nucleotides present within the given sequence
checklist =[]
for x in nucleotides:
    if x not in checklist:
        checklist.append(x)
print(checklist)
"""

print ("Raw Counts")
print ("A: ", num_a)
print ("C: ", num_c)
print ("G: ", num_g)
print ("T: ", num_t)
print ("N: ", num_n)
print ("Total: ", total_nuc)

## Part 3
### TODO Print out the frequencies for each nucleotide in alphabetical order
#Dividing the individual nucleotides (A, C, G, T) by the total number of nucleotides, turning the float into a decimal and then rounding to the hundredth
"""
total_nuc_freq = num_a + num_c + num_g + num_t
"""

freq_a = round(Decimal(num_a/total_nuc), 3)
freq_c = round(Decimal(num_c/total_nuc), 3)
freq_g = round(Decimal(num_g/total_nuc), 3)
freq_t = round(Decimal(num_t/total_nuc), 3)

print ("\nNucleotide Frequencies")
print ("A: ", freq_a)
print ("C: ", freq_c)
print ("G: ", freq_g)
print ("T: ", freq_t)

## Part 5
### TODO Use overlapping windows to count the dinucleotides. See the assignment for more information on overlapping wiindows.

#Intializing a dictionary and list
dinucleotide_dictionary = {}
dinuc_list = []

#Populating my list with "recognized" dinucleotides consisting of (A, C, G, T) and not (N, R, W, Y) 
for nuc1 in ['A', 'C', 'G', 'T']:
    for nuc2 in ['A', 'C', 'G', 'T']:
        dinucleotide = nuc1 + nuc2
        dinuc_list.append(dinucleotide)

#Iterating through the given nucleotide sequence and adding any dinucleotide to my dictionary from the sequence given
for i in range(0, len(nucleotides)-1):
    next_dinuc = nucleotides[i] + nucleotides[i+1]
    if next_dinuc in dinucleotide_dictionary:
        dinucleotide_dictionary[next_dinuc] += 1
    else:
        dinucleotide_dictionary[next_dinuc] = 1

#Intializing a list and filling it with all the keys from my dictionary
all_dinuc_list = list(dinucleotide_dictionary.keys())

#Iterating through all_dinuc_list
for k in all_dinuc_list:
#I remove all dinucleotides from my dictionary if they do not exist in my "recognized" dinucleotide list
    if k not in dinuc_list:
        dinucleotide_dictionary.pop(k)

#Summing all values of the keys within my dinucleotide_dictionary
total_dinuc = sum(dinucleotide_dictionary.values())

#Calculating all the frequencies for the keys within dinucleotide_dictionary
for value in dinucleotide_dictionary:
    dinucleotide_dictionary[value] = round(Decimal(dinucleotide_dictionary[value]/total_dinuc), 3)

#Sorting my dinucleotide_dictionary in alpehbetical order
ordered = collections.OrderedDict(sorted(dinucleotide_dictionary.items()))

### TODO Print the frequencies for each of dinucleotides in alphabetical order
print("\nDinucleotide Frequencies")
for key in ordered:
    print(key, ':', ordered[key])
