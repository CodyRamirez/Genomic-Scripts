Question 1:
Command to run genotype_analysis.py on the sperm VCF file:
python3 genotype_analysis.py --vcf /home/assignments/assignment12/sperm_genotype_calls.vcf

Refer to Table 1) for:
~ The mean and variance in variant coverage for each dataset
~ Variants that have unknown coverage for each dataset

Table 1)
Dataset_Name	Mean	Variance	Unknown	Allelic_Dropout	Allelic_Dropout_Rate
Blood Standard	22.15	241.46631	2517	'TRUTH'	'TRUTH'
Live Standard	70.33	1833.8289	599		1214	0.08
Dead Standard	95.27	13017.401	47905	6143	0.4
Live MALBAC1	80.42	11001.722	20600	6476	0.42
Live MALBAC2	56.77	6090.1249	41589	6141	0.4
Live MALBAC3	49.59	5169.6695	18204	5645	0.37
Dead MALBAC1	78.72	2497.4378	689		1264	0.08
Dead MALBAC2	79.92	8716.1191	20465	6915	0.45
Dead MALBAC3	59.54	6483.3149	13842	5287	0.34

The mean coverage varies amoung each dataset, the mean variant coverage ranges from blood standard mean of 22.15 and dead standard mean of 95.27. The variance varies amoung each dataset, the variance ranges from blood standard variance of 241 and dead standard variance of 13,017. This is not what I expected, I had expected the mean coverage to be higher and the variance lower for the MALBAC datasets.
-
Question 2:
 Refer to Table 1) column labeled Allelic_Dropout_Rate for the rate of allelic dropout for each of the non-blood datasets.
 The numerators are located in Table 1) column labeled Allelic_Dropout.
 The number of heterozygous sites identified in the blood standard dataset (the denominator): 15376
-
Question 3
Refer to Table 2) for:
~ The number of variants in the list of “dead” variants
~ The number of variants in the list of 'live' variants

Table 2)
Variant	Count	TransI	TransV	Transitions/Transversions
'Live'	1321	797		524		1.52
'Dead'	501		288		213		1.35

-
Question 4a
Refer to Table 2) for:
~ The number of transitions in the list of “dead” variants
~ The number of transversions in the list of “dead” variants
~ The transition/transversion ratio in the list of “dead” variants

-
Question 4b
Refer to Table 2) for:
~ The number of transitions in the list of “live” variants
~ The number of transversions in the list of “live” variants
~ The transition/transversion ratio in the list of “live” variants

-
Question 4c
Yes; count, transitions and transversions are different between the two sets. It could mean that more viable genetic information to sequence existed for the 'live' dataset instead of the 'dead' dataset even though the transition/transversion ratio is approximately the same.
-
Question 5
I can conclude that transversions occur more often within 'dead' versus 'live' sperm, which would make sense considering that transversions are more likely to cause nonsynonymous mutations rather than synonymous.
-
Question 6
I can conclude that genotype calls from small pools of data suffer more from PCR bias compared to bulk considering the variance and that artifacts should be considered.
-
Extra credit
{Suggest three other ideas for characterizing the mutational differences between live and dead sperm cells.}
-
Comments:
{Things that went wrong or you can not figure out}
-
Suggestions:
{What programming and/or genomics topics should the TAs cover in the next class that would have made this assignment go smoother?}