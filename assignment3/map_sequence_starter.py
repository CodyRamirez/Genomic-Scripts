#!/usr/bin/env python3
"""
Align sequencing reads to a set of genes.

Usage: map_sequence_starter.py <cDNA_fasta> <seq_reads_file>

Arguments:
    <cDNA_fasta>     = Path to cDNA fasta file.  (Required.)
    <seq_reads_file> = Path to sequencing reads file.  (Required.)
"""

#Import modules needed for program
import sys

def reverse_complement(sequence):
    """This function takes the reverse complement of a sequence"""

    #define complement dictionary
    complement_dictionary = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N' }

    #Reversing the string by using extended slice syntax which works by doing "[begin:end:step]", by leaving begin and end out and specifying a step of -1, it reverses a strings order.
    sequence = sequence[::-1]
    
    #Taking the complement by turning the sequence into a list of bases so I can turn a single base at a time into it's complement base using a list comprehension
    bases = list(sequence)
    bases = [complement_dictionary[base] for base in bases]#This takes a base from the bases list and switches it for its complement base

    #return the reverse complement by joining all of the bases from the list bases
    return ''.join(bases)


def read_cDNA_file_to_dict(filename):
    """This function reads a cDNA file into a dictionary."""
    
    #initialize dictionary
    cDNA_dictionary = {}

    #open file
    with open(cDNA_file) as f:
    
        #loop through file line by line
        for line in f:

            #remove newline
            line = line.rstrip()
        
            #get gene name
            if line.startswith(">"):#If the line starts with the character ">" then,
                gene_name = line.split("|")[1]#I separate the line by the character "|" and assign index 1 to gene_name
        
            #read in sequence in uppercase
            if not line.startswith(">"):#If the line does not start with the character ">" then,
                line = line.upper()#I make all of the characters within the line uppercase

            #put name and sequence in dictionary
            cDNA_dictionary[gene_name] = line#I assign the gene_name as the key and the line (sequence) as the value

    #return dictionary    
    return cDNA_dictionary


def create_twenty_five_mer_dict(cDNA_dict):
    #This function creates a dictionary of 25mer sequences from the cDNA sequences.  The dictionary keys are the 25mer sequences and the values are the corresponding gene names.   

    #initialize dictionary
    twenty_five_mer_dictionary = {}

    #Looping through the keys of the cDNA dictionary which are the gene_names with corresponding sequence as their value
    for key in cDNA_dict:

        #get sequence that corresponds to the key
        sequence = cDNA_dict[key]

        #move through sequence, grabbing 25 mers and create dictionary
        for i in range(0, len(sequence)-24):#Navigating through the sequence an index (base) at a time
                current_25mer = sequence[i:i+25]#Grabbing a 25mer by using the current index to create a 25 base range and assigning that to current_25mer
                if current_25mer in twenty_five_mer_dictionary:#If the current_25mer exist within my 25mer-dictionary then,
                    twenty_five_mer_dictionary[current_25mer] = key#I assign the 25mer as the key and the gene_name its value
                else:
                    twenty_five_mer_dictionary[current_25mer] = key#Else, I assign the 25 mer as the key and the gene_name its value

    #return dictionary
    return twenty_five_mer_dictionary

##############main loop###################

#parse command line
#check that the correct number of arguments were given
if (len(sys.argv) != 3):
    sys.exit(__doc__)

cDNA_file = sys.argv[1]
seq_reads_file = sys.argv[2]

#read in cDNAs
cDNA_dict = read_cDNA_file_to_dict(cDNA_file)
print("Read in the cDNAs")

#read in 25mers
twenty_five_mer_dict = create_twenty_five_mer_dict(cDNA_dict)
print("Created dictionary")


count_dict = {}#Creating a count dictionary
for key in cDNA_dict:#Loops through every key in cDNA_dict (keys are the gene_names)
    count_dict[key] = 0#Adds every key from cDNA_dict as a key to count_dict with the value intial value of 0


fh = open(seq_reads_file, "r")#Opening the file to read

for line in fh:#Looping through each line

    line = line.rstrip()#Returns a copy of line in which all whitespace characters have been stripped from the end of the string
    line = line.upper() #Returns a copy of line in which all case-based characters have been uppercased
    if line in twenty_five_mer_dict.keys():#If the current line is found within the twenty_five_mer_dict key list then,
        if twenty_five_mer_dict[line] in count_dict.keys():#If the current lines corresponging gene_name also exist within count_dict then,
            count_dict[twenty_five_mer_dict[line]] = count_dict[twenty_five_mer_dict[line]] + 1#Increment that gene_name by one within count_dict

    else: 
        rc_seq = reverse_complement(line)#Assigns the reverse complement of the current line to rc_seq
        if rc_seq in twenty_five_mer_dict.keys():#If the rc_seq exist within the twenty_five_mer_dict then,
            if twenty_five_mer_dict[rc_seq] in count_dict.keys():#If the current rc_seq gene_name exist within count_dict then,
                count_dict[twenty_five_mer_dict[rc_seq]] = count_dict[twenty_five_mer_dict[rc_seq]] + 1#Increment that gene_name by one within count_dict

fh.close()#Close the file of seq_reads_file that was being read

print("Name\tReads\tReads per BP")#Printing the tab delimited line
for name, count in sorted(count_dict.items()):#Looping through the sorted keys and values of count_dict
    print(name + "\t" + str(count) + "\t" + str(float(count)/len(cDNA_dict[name])))#Printing the gene_name, count of gene_name and the number of gene_name counts divided by the length of the sequence to compisate for the size
