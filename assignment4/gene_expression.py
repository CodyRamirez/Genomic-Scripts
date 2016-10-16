#!/usr/bin/env python3
"""
This script calculates the expression of genes given raw counts of RNA-seq as input from a table where genes are the rows and samples are the columns
Usage: python3 gene_expression.py raw_counts.txt
"""


import sys
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


#TODO: Print out the doc string and exit if the number of input parameters is not correct
#sys.arg is a list containing 2 elements: the script name and a command line argument
#Checks that the input paramenters are correct. If not, prints an error statement, documentation and exit.
if (len(sys.argv)!=2):
    sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)


#######################
##Part 0 -- Functions##
#######################


# CPM counts per million
# Convert raw counts to counts per million (cpm)
# Raw count x[i,j] (from sample i, gene j)
# Total counts N[i] (from sample i)
# cpm[i,j] = (10^6)*x[i,j]/N[i]
def counts_per_million(dictionary, list_of_samples):
    # Initialize output dictionary cpm_dict
    cpm_dict = {}
    # Find N, the list of each sample's library size
    N = library_sizes(dictionary, list_of_samples)
    # Calculate cpm one gene at a time using list comprehension
    for k,v in dictionary.items():
        # Note: do NOT use 10^6, which equals 12
        # Use either 10**6 or 1000000 for 1 million
        # k is the key (gene name)
        # v is the value (raw counts of gene k)
        cpm_dict[k] = [(10**6)*x/n for x,n in zip(v,N)]
    return(cpm_dict)
# End of CPM function #


# Library size
# Calculate library size of each sample (e.g. sum of RNA-seq counts)
def library_sizes(dictionary, list_of_samples):
    # Get the total number of samples
    num_samples = len(list_of_samples)
    # Initialize N, a list to hold the total counts from each sample
    N = []
    # For loop to iterate over each sample 
    for i in range(num_samples):
        # Append a new float zero value for each sample (goes to index i)
        N.append(0.0)
        # For loop to iterate over each value in our dictionary
        for v in dictionary.values():
            # Get the count from the i index of this gene and add it to the total for sample i
            N[i] += v[i]
    # Return the list containing each sample's library size
    return(N)
# End of lbirary sizes function #


# TODO: Fill in the code for the translate_dictionary function.
# Translate dictionary
# Function to translate a dictionary from {gene:[list of counts by sample]} to {sample:[list of counts by gene]}
# The new dictionary will have one key for each sample and the value of each key will be a list of counts associated with that sample
# This function is used in another function called upper_quartile_norm()
# Each comment below should correspond to one line of code
def translate_dictionary(dictionary, list_of_samples):
    # Initialize a new dictionary called translated_dictionary with empty curly brackets
    translated_dictionary = {}
    # Sort the keys of the input dictionary and save it in a list called 'genes'
    genes = sorted(list(dictionary))
    # Use a for loop to iterate over each sample (using a numbered index, not sample names)
    for i in range(len(list_of_samples)):
        # Get the current sample name from list_of_samples
        current_sample_name = list_of_samples[i]
        # Initialize a new key in translated_dictionary using the current sample name. Let the value be an empty list.
        translated_dictionary[current_sample_name] = []
        # Use a for loop to iterate over each gene in your list of genes
        for j in genes:
            # Find the RNA-seq count associated with this gene for this sample
            rna_seq_count = dictionary[j][i]
            # Append the count to the list of counts for this sample in translated_dictionary
            translated_dictionary[current_sample_name].append(rna_seq_count)
    # Return translated_dictionary
    return(translated_dictionary)
# End of translate dictionary function #


# TODO: Add comments to the code in the upper_quartile_norm function to explain what each line does
# Upper quartile normalization
# Compute the upper quartile normalization of raw counts
# Raw count x[i,j] (from sample i, gene j)
# D[i] value corresponding to 75th percentile of raw counts (from sample i)
# Mean of D, or the mean of all upper quartile means
# Upper quartile normalized count for sample i and gene j is (Mean of D)*x[i,j]/D[i]
def upper_quartile_norm(dictionary, list_of_samples):
    # Initialize a new dictionary called upper_quartile_norm_dictionary with empty curly brackets 
    upper_quartile_norm_dictionary = {}
    # Initialize a variable num_samples that is equal to the number of samples in the list list_of_samples
    num_samples = len(list_of_samples)
    # Initialize a new dictionary called translated_dict_for_D equal to the transposition of the current dictionary
    translated_dict_for_D = translate_dictionary(dictionary, list_of_samples)
    # Initialize a new list called D with empty brackets
    D = []
    # Use a for loop to iterate through the num_samples variable 
    for i in range(num_samples):
        # Assign sample_name as the current sample name indexing through list_of_samples 
        sample_name = list_of_samples[i]
        # Appending the 75th percentile of raw counts from the current sample_name to list D 
        D.append(np.percentile(translated_dict_for_D[sample_name],75))
    # Assigning meanD to the mean of the list D which holds the values of the 75th percentile of raw counts from all samples 
    meanD = np.mean(D)
    # Use a for loop to iterate through the corresponging keys and values of the dictionary 
    for k,v in dictionary.items():
        # Assigning the current gene (k) as the key for upper_quartile_norm_dictionary and the upper quartile normalized count for sample i and gene j as the value
        upper_quartile_norm_dictionary[k] = [meanD*x/d for x,d in zip(v,D)]
    # Returns upper_quartile_norm_dictionary
    return(upper_quartile_norm_dictionary)
# End of upper quartile normalization function #


# TODO: Write a function to calculate Fisher's Linear Discriminant (add comments, too!) for all genes in your count dictionary. Call your function fishers_linear_discriminant. Remember that functions should be able to work using different data sets, so make sure that your function would work with data from a different experiment with different numbers of before and after samples. 
# Return the top ten genes according to the highest FLD values
# Input to this function should be a count dictionary and two lists letting the function know which index values go with each group
# You will need to iterate over each gene in the dictionary and calculate the FLD of each gene
# FLD of a gene = ((m1-m2)^2)/((s1)^2 + (s2)^2)
# m1 = mean of the first group
# m2 = mean of the second group
# s1 = standard deviation of the first group
# s2 = standard deviation of the second group
def fishers_linear_discriminant(dictionary, list1, list2):
    # Initializing a new empty dictionary called FLD_dictionary to hold gene:FLD
    FLD_dictionary = {}
    # Use a for loop to iterate through each gene in the dictionary passed in
    for k,v in dictionary.items():
        # Initialize two list: counts1 and counts2
        counts1 = []
        counts2 = []
        # Use a for loop to iterate through list1
        for i in list1:
            # Appending the current value of the dictionary with the ith index (based on list1) to counts1 list
            counts1.append(v[i])
        # Use a for loop to iterate through list2
        for j in list2:
            # Appending the current value of the dictionary with the jth index (based on list2) to counts2 list 
            counts2.append(v[j])
        # Assigning m1 to the mean of the first group and m2 to the mean of the second group
        m1 = np.mean(counts1)
        m2 = np.mean(counts2)
        # Assigning s1 to the standard deviation of the first group and s2 to the standard deviation of the second group
        s1 = np.std(counts1)
        s2 = np.std(counts2)
        # Calculateing the FLD of the current gene
        current_FLD = ((m1-m2)**2)/((s1)**2 + (s2)**2)
        # Adding gene as the key and it's corresponging FLD as it's value 
        FLD_dictionary[k] = current_FLD
    # Taking the top ten genes according to the highest FLD values and putting them in the top_ten list    
    top_ten = sorted(FLD_dictionary, key=FLD_dictionary.get, reverse=True)[:10]
    # Use a for loop to iterate through each gene in the top_ten list
    for gene in top_ten:
        # Prints out the top_ten genes and corresponding FLD values
        print(gene + ':' + str(FLD_dictionary[gene]))
# End of fishers linear discriminant function #


####################
##End of functions##
####################


# These first lines of code will get the data imported and in the right format for the rest of the homework
# You need to do the rest of the work starting from Part 1 -- Data filtering

#open the data file
data_file = open(sys.argv[1],"r")
#first line of data file contains sample names -- store in a list
sample_list = data_file.readline().strip().split()[1:]
#initialize raw count dictionary
raw_counts_dict = {}
#add each gene and expression values to dictionary
for line in data_file:
    line_list = line.strip().split() #split line into list at whitespace, strip removes leading and trailing whitespace
    gene = line_list[0] #name of gene is the first thing in the list
    expression_values = [int(float(v)) for v in line_list[1:]] #values of gene expression follow the gene name
    raw_counts_dict[gene] = expression_values #add keys and values to the dictionary {gene:expression}
#close the data file
data_file.close()


############################
##Part 1 -- Data filtering##
############################


# Filter out genes with zero expression in all samples
# Initializing a new dictionary called genes_pass_first_filter which will hold all genes with their values
genes_pass_first_filter = {}
# Initializing a new dictionary called genes_fail_first_filter which will hold all genes with their values 
genes_fail_first_filter = {}
# Initializing genes_pass_first_filter_count and genes_fail_first_filter_count to 0
genes_pass_first_filter_count = 0
genes_fail_first_filter_count = 0

# Use a for loop to iterate through the key:value (Note: the value in the dictionary is a list of values) pairs of raw_counts_dict
for k,v in raw_counts_dict.items():
    # Use a IF statement to ensure that the values do not sum 0
    if sum(v) != 0:
        # Assign the current gene and it's values to genes_pass_first_filter
        genes_pass_first_filter[k] = v
        # Increment genes_pass_first_filter by 1
        genes_pass_first_filter_count += 1
    else:
        # Use a ELSE statement to ensure that the current gene and it's values that sum to 0 are assigned to genes_fail_first_filter
        genes_fail_first_filter[k] = v
        # Increment genes_fail_first_filter by 1
        genes_fail_first_filter_count += 1
# Prints out the number of genes that passed the first filter
print("There are " + str(genes_pass_first_filter_count) + " genes that pass the first filter!")


# Filter out genes with 20 or more samples with cpm < 1
# Initializing a new dictioanry called genes_pass_second_filter which will hold all gene with 20 or more samples with cpm > 1
genes_pass_second_filter = {}
# Initializing a new dictionary called genes_fail_second_filter which will hold all gene with 20 or more samples with cpm < 1
genes_fail_second_filter = {}
# Initializing genes_pass_second_filter_count and genes_fail_second_filter_count to 0
genes_pass_second_filter_count = 0
genes_fail_second_filter_count = 0
# Assigning cpm_dict to counts_per_million of genes_pass_first_filter
cpm_dict = counts_per_million(genes_pass_first_filter, sample_list)

# Use a for loop to iterate through genes_pass_first_filter
for k,v in genes_pass_first_filter.items():
    # Initializing cpm_fail_count to 0
    cpm_fail_count = 0
    # Use a for loop to iterate through the cpm_values in cpm_dict
    for cpm_value in cpm_dict[k]:
        # Use a IF statement to determine if the current cpm_value is less than one for the current gene (k)
        if cpm_value < 1:
            # Increment cpm_fail_count by 1
            cpm_fail_count += 1
    # Use a IF statement to determine if the current cpm_fail_count is greater or equal to 20 than
    if cpm_fail_count >= 20:
        # Assign the current gene and it's corresponding values to genes_fail_second_filter
        genes_fail_second_filter[k] = v
        # Increment genes_fail_second_filter_count by 1
        genes_fail_second_filter_count += 1
    else:
        # Assign the current gene and it's corresponding values to genes_pass_second_filter
        genes_pass_second_filter[k] = v
        # Increment genes_pass_second_filter_count by 1
        genes_pass_second_filter_count += 1
# Print out the number of genes_pass_second_filter_count
print("There are " + str(genes_pass_second_filter_count) + " genes that pass the second filter!")


################################
##Part 2 -- Data visualization##
################################


# Plot library sizes (save file as library_size.png)

# Assigning the library size of genes_pass_second_filter to library_size
library_size = library_sizes(genes_pass_second_filter, sample_list)
# Print the range of the library size
print('The range of the library size after filtering is: (min = ' + str(min(library_size)) + ', max = ' + str(max(library_size)) +')')

fig = plt.figure()# Start a new figure
plt.bar(left = range(40), height = library_size, align = 'center')# Making a bar plot
plt.xticks(range(len(sample_list)), sample_list, size = 'small', rotation = 'vertical')# Setting the x-limits of the current tick locations and labels
plt.xlabel('Samples')# Assigning the x-axis label
plt.ylabel('Total million counts')# Assigning the y-axis label
plt.title('Total million counts of samples')# Assigning the title
plt.tight_layout# Automatically adjust plot parameters to give specified padding
fig.savefig('library_size.png')# Saving the current bar plot as library_siz.png


################################
##Part 3 -- Data normalization##
################################


# Normalize count data left after filtering steps
normalized_counts_dictionary = upper_quartile_norm(genes_pass_second_filter, sample_list)

# Plot normalized library sizes (save file as library_size_normalzied.png)
# Assigning the library size of normalized_counts_dictionary to normalized_library_size
normalized_library_size = library_sizes(normalized_counts_dictionary, sample_list)
# Print the range of normalized_library_size
print('The range of the library size after filtering is: (min = ' + str(min(normalized_library_size)) + ', max = ' + str(max(normalized_library_size)) +')')

fig = plt.figure()# Start a new figure
plt.bar(left = range(40), height = normalized_library_size, align = 'center')# Making a bar plot
plt.xticks(range(len(sample_list)), sample_list, size = 'small', rotation = 'vertical')# Setting the x-limits of the current tick locations and labels
plt.xlabel('Samples')# Assigning the x-axis label
plt.ylabel('Total million counts')# Assigning the y-axis label
plt.title('Total million counts of samples normalized')# Assigning the title
plt.tight_layout# Automatically adjust plot parameters to give specified padding
fig.savefig('library_size_normalzied.png')# Saving the current bar plot as library_size_normalized.png


##############################
##Part 4 -- Data exploration##
##############################


# Calculate Fisher's Linear Discriminant for each gene basd on the normalized count data
fishers_linear_discriminant(normalized_counts_dictionary, range(20), range(20,40))

# Pick a gene to explore further and plot the mean expression level for the before and after groups (save file as mean_expression.png)
# Initializing top_gene_normalized_counts_dictionary
top_gene_normalized_counts_dictionary = {}
# Assiging gene GRB14 and it's normalized counts to top_gene_normalized_couts_dictionary
top_gene_normalized_counts_dictionary['GRB14'] = normalized_counts_dictionary['GRB14']
# Initializing a list groups with Before and After
groups = ['Before', 'After']
# Initializing two list: group1 and group2
group1 = []
group2 = []
# Assigning the first 20 samples (the before samples) to group1
group1 = top_gene_normalized_counts_dictionary['GRB14'][0:20]
# Assigning the last 20 samples (the after samples) to group2
group2 = top_gene_normalized_counts_dictionary['GRB14'][20:40]
# Calculating the mean values for group1 and group2
group1_mean = np.mean(group1)
group2_mean = np.mean(group2)
# Assigning the group mean values to groups_mean
group_means = (group1_mean, group2_mean)
# Calculating the standard deviation values for group1 and group2
group1_std = np.std(group1)
group2_std = np.std(group2)
# Calculating the standard error of the mean (SEM) for group1 and group2
group1_sem = (group1_std)/np.sqrt(20)
group2_sem = (group2_std)/np.sqrt(20)
# Assigning the group SEM values to group_sems
group_sems = (group1_sem, group2_sem)

fig = plt.figure()# Start a new figure
plt.bar(left = range(2), height = group_means, width = 0.5, yerr = group_sems, align = 'center')# Making a bar plot
plt.xticks(range(2), groups, size = 'small', rotation = 'vertical')# Setting the x-limits of the current tick locations and labels
plt.xlabel('Groups')# Assigning the x-axis label
plt.ylabel('Mean Expression Level')# Assigning the y-axis label
plt.title('Mean Expression Level for Gene: GRB14')# Assigning the title
plt.tight_layout# Automatically adjust plot parameters to give specified padding
fig.savefig('mean_expression.png')# Saving the current bar plot as mean_expression.png

"""
######################################
##Extra credit -- Going the distance##
######################################

# Calculate Euclidean distance matrix of samples, output the most/least related samples

# Plot a dendrogram using Euclidean distance (save file as dendrogram.png)
"""