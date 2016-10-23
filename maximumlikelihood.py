#!/usr/bin/env python3
"""
Assumes that the quantitative phenotype is normally distributed within each of the 3 unknown disease genotype groups.
Assumes that all parents are statistically independent (unrelated).
It cannot use the disease genotype (which are unknown) and the marker genotype are irrelevant.

INPUT: This script uses as input data the quantiative phenotype (QTp) on the parents only (ignoring the children's data).

OUTPUT: The program calculates and outputs the maximum likelihood estimator of 5 parameters:
	mudaa: Mean phenotype for all subjects with disease genotype A/A
	mudab: Mean phenotype for all subjects with disease genotype A/B
	mudbb: Mean phenotype for all subjects with disease genotype B/B
	de: Common variance (standard deviation) of the quantiative phenotype within each of the 3 disease genotype groups
	q: Allele frequency of (one of) the disease gene alleles

Usage: python3 maximumlikelihood.py <list of QTp values> <number of grid steps>
"""


import os
import sys
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import timeit


########## STARTING TIMER ##########
start_time = timeit.default_timer()
####################################


# Check if an argument was passed to the python script
# If not exit exit and print the documentation
if (len(sys.argv) != 3):
	sys.exit("ERROR: incorrect number of arguments.\n" + __doc__)

# Assign input variable to QTP_file
QTp_file = sys.argv[1]
grid_steps = int(sys.argv[2])

# Creating dynamic output filenames by grabbing the filename and parsing out the basename
basename = os.path.basename(sys.argv[1])
# Remove the extension
basename_no_extension = os.path.splitext(basename)[0]

# Initializing QTp_list to hold all QTp values
QTp_list = []
# Initializing alt_current_max to hold the current maximum likelihood
alt_current_max = 0
# Initializing alt_parameters_list to hold all parameters that created the maximum likelihood probability
alt_parameters_list = []
# Initializing population_likelihood to hold maximumlikelihood value
population_likelihood = 0

#Opening the QTp_file as f
with open(QTp_file) as f:
	for line in f:# Using a for loop to iterate through every line in QTp_file
		line_list = line.strip().split()# Split line into list at whitespace, strip removes leading and trailing whitespace
		QTp_list.append(float(line_list[0]))# Appending the current QTp value to the QTp_list


mu_stepsize = (max(QTp_list)-min(QTp_list))/grid_steps
# Testing the maximumlikelihood for the alternative hypothesis
for mudaa in np.arange(min(QTp_list), max(QTp_list) + mu_stepsize, mu_stepsize):
	for mudab in np.arange(min(QTp_list), max(QTp_list) + mu_stepsize, mu_stepsize):
		for mudbb in np.arange(min(QTp_list), max(QTp_list) + mu_stepsize, mu_stepsize):
			for de in np.arange(0.001, np.std(QTp_list), (np.std(QTp_list))/grid_steps):
				for q in np.arange(0, 0.6, 0.1):

					mudaa = round(mudaa, 2)
					mudab = round(mudab, 2)
					mudbb = round(mudbb, 2)
					q = round(q, 2)

					current_likelihood = 0
					natural_log_of_current_likelihood = 0
					population_likelihood = 0
					
					for QTp in QTp_list:
						current_likelihood = ((1-q)**2)*(1/np.sqrt(2*(np.pi)*(de)))*(np.exp((-0.5)*(((QTp-mudaa)/de)**2))) + (2*q)*(1-q)*(1/np.sqrt(2*(np.pi)*(de)))*(np.exp((-0.5)*(((QTp-mudab)/de)**2))) + (q**2)*(1/np.sqrt(2*(np.pi)*(de)))*(np.exp((-0.5)*(((QTp-mudbb)/de)**2)))
						natural_log_of_current_likelihood = np.log1p(current_likelihood)
						population_likelihood += natural_log_of_current_likelihood

						if population_likelihood > alt_current_max:
							alt_current_max = population_likelihood
							alt_parameters_list = [mudaa, mudab, mudbb, de, q]
							#print('The current max is: ' + str(alt_current_max) + ' using these parameters: ' + str(alt_parameters_list))
						#else:
							#print('Population likelihood:' + str(population_likelihood) + ' < current max:' + str(alt_current_max) + ' max parameters: ' + str(alt_parameters_list))						
print('minimum QTp value: ' + str(min(QTp_list)))
print('maximum QTp value: ' + str(max(QTp_list)))

print('\nThe maximum likelihood for the alternative hypothesis is: ' + str(alt_current_max) )
print('These parameters were used to calculate the maximum likelihood:')
print('mudaa: ' + str(alt_parameters_list[0]))
print('mudab: ' + str(alt_parameters_list[1]))
print('mudbb: ' + str(alt_parameters_list[2]))
print('de: ' + str(alt_parameters_list[3]))
print('q: ' + str(alt_parameters_list[4]))



null_current_max = 0
null_parameters_list = []
# Testing for the null hypothesis
for mudaa in np.arange(min(QTp_list), max(QTp_list) + mu_stepsize, mu_stepsize):
	for de in np.arange(0.001, np.std(QTp_list) + 0.001, (np.std(QTp_list))/grid_steps):
		for q in np.arange(0, 1.1, 0.1):

			mudaa = round(mudaa, 2)
			mudab = round(mudab, 2)
			mudbb = round(mudbb, 2)
			q = round(q, 2)


			current_likelihood = 0
			natural_log_of_current_likelihood = 0
			population_likelihood = 0
			
			for QTp in QTp_list:
				current_likelihood = ((1-q)**2)*(1/np.sqrt(2*(np.pi)*(de)))*(np.exp((-0.5)*(((QTp-mudaa)/de)**2))) + (2*q)*(1-q)*(1/np.sqrt(2*(np.pi)*(de)))*(np.exp((-0.5)*(((QTp-mudaa)/de)**2))) + (q**2)*(1/np.sqrt(2*(np.pi)*(de)))*(np.exp((-0.5)*(((QTp-mudaa)/de)**2)))
				natural_log_of_current_likelihood = np.log1p(current_likelihood)
				population_likelihood += natural_log_of_current_likelihood

				if population_likelihood > null_current_max:
					null_current_max = population_likelihood
					null_parameters_list = [mudaa, mudaa, mudaa, de, q]
					#print('The current max is: ' + str(null_current_max) + ' using these parameters: ' + str(null_parameters_list))
				#else:
					#print('Population likelihood: ' + str(population_likelihood) + ' is not greater than the current max: ' + str(null_current_max))
print('\nThe maximum likelihood for the null hypothesis is: ' + str(null_current_max) )
print('These parameters were used to calculate the maximum likelihood:')
print('mudaa: ' + str(null_parameters_list[0]))
print('mudab: ' + str(null_parameters_list[1]))
print('mudbb: ' + str(null_parameters_list[2]))
print('de: ' + str(null_parameters_list[3]))
print('q: ' + str(null_parameters_list[4]))



chi_square_value = -2*null_current_max - (-2*alt_current_max)
p_value = stats.chisqprob(chi_square_value, 3)
if p_value <= 0.05:
	print('\nReject the Null Hypothesis because the p-value = ' + str(p_value))
else:
	print('\nFailed to reject the null hypothesis because the p-value = ' + str(p_value))



min_QTp = min(QTp_list)
max_QTp = max(QTp_list)
alt_mudaa = alt_parameters_list[0]
alt_mudab = alt_parameters_list[1]
alt_mudbb = alt_parameters_list[2]
alt_de = alt_parameters_list[3]

histogram_plot_name = basename_no_extension + '_histogram_plot.png'# Create histogram plot name
fig = plt.figure()# Start a new figure
plt.hist(QTp_list, bins = 100, histtype = 'bar', align = 'mid')# Plotting QTp values using 100 bins with bars aligned to the middle
plt.xlabel('Phenotype values')# X-axis label
plt.ylabel('Counts')# Y-axis label
plt.title('Quantitative Phenotype (QTp) Distribution')# Title of histogram
fig.savefig(histogram_plot_name)# Saving figure as QTp_plot

density_QTp_plot_name = basename_no_extension + '_density_plot.png'# Create density plot name
fig = plt.figure()# Start a new figure
data = QTp_list
density = stats.kde.gaussian_kde(data)
x = np.arange(min_QTp-10, max_QTp+10, 1)
density.covariance_factor = lambda : .25
density._compute_covariance()
plt.plot(x, density(x))
plt.xlabel('Phenotype values')# X-axis label
plt.ylabel('Density')# Y-axis label
plt.title('Quantitative Phenotype (QTp) Distribution')# Title of histogram
fig.savefig(density_QTp_plot_name)# Saving figure as QTp_plot

normal_distribution_plot_name = basename_no_extension + '_normal_distribution_plot.png'# Create normal distribution plot name
fig = plt.figure()
x1 = np.linspace(alt_mudaa - alt_de*10, alt_mudaa + alt_de*10, 100)
x2 = np.linspace(alt_mudab - alt_de*10, alt_mudab + alt_de*10, 100)
x3 = np.linspace(alt_mudbb - alt_de*10, alt_mudbb + alt_de*10, 100)
y1 = stats.norm.pdf(x1, loc = alt_mudaa, scale = alt_de)
y2 = stats.norm.pdf(x2, loc = alt_mudab, scale = alt_de)
y3 = stats.norm.pdf(x3, loc = alt_mudbb, scale = alt_de)
plt.plot(x1, y1)
plt.plot(x2, y2)
plt.plot(x3, y3)
fig.savefig(normal_distribution_plot_name)


########## ENDING TIMER ##########
stop_time = timeit.default_timer()
##################################
print('\nTime to run program in minutes: ' + str(((stop_time - start_time)/60)))




