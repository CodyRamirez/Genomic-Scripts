#!/usr/bin/env python3

"""make_seq.py prints a random sequence given a sequence length and nucleotide frequencies.  The random sequence will have the same nucleotide frequencies as the input nucleotide frequencies.

Usage: python3 make_seq.py <sequence_length> <a_freq> <c_freq> <g_freq> <t_freq>

<sequence_length> = Length of sequence to generate
<a_freq> = frequency of As
<c_freq> = frequency of Cs
<g_freq> = frequency of Gs
<t_freq> = frequency of Ts
"""

# Import modules 
import sys
import random
from decimal import Decimal 

# sys.arg is a list containing 6 elements: the script name and 5 command line arguments
# Check that all 5 command line arguments were given. If not, print the documentation and exit.
if (len(sys.argv) != 6):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__) 

# Save the input arguments as variables
# By default, the command line arguments are saved as strings. Convert them to numeric types.
length = int(sys.argv[1])
a_freq = float(sys.argv[2])
c_freq = float(sys.argv[3])
g_freq = float(sys.argv[4])
t_freq = float(sys.argv[5])

# Check that frequencies add to 1. If not, exit the program 
if (abs(a_freq + t_freq + c_freq + g_freq - 1) > 1e-4):
	sys.exit("ERROR: Nucleotide frequencies do not add up to 1!")

## Part 4
### TODO Generate a random nucleotide sequence

# Initialize an empty string that nucleotides can be appended to
seq = ""

if a_freq == 0:
    num_a = 0
else:
    num_a = round(round(Decimal(length*a_freq), 1))

if c_freq == 0:
    num_c = 0
else:
    num_c = round(round(Decimal(length*c_freq), 1))

if g_freq == 0:
    num_g = 0
else:
    num_g = round(round(Decimal(length*g_freq), 1))

if length != (num_a + num_c + num_g):
    num_t = round(length-(num_a + num_c + num_g))
else:
    num_t = 0

# Create a for loop that will be repeated <length> times
#for i in range(0, length):
while (length != len(seq)):
#     Generate a random decimal
    rand = random.random()

    next_nuc = ""

#     Use if/else if/else logic to determine which nucleotide to add
    if rand <= 0.25 and num_a != 0:
        num_a = num_a - 1
        next_nuc = "A"
    elif rand > 0.25 and rand < 0.5 and num_c != 0:
        num_c = num_c - 1
        next_nuc = "C"
    elif rand >= 0.5 and rand < 0.75 and num_g != 0:
        num_g = num_g - 1
        next_nuc = "G"
    elif rand >= 0.75 and num_t != 0:
        num_t = num_t - 1
        next_nuc = "T"
#     Append the nucleotide to the nucleotide sequence
    seq += next_nuc
    next_nuc = ""

# Print the full nucleotide sequence
print (seq)
