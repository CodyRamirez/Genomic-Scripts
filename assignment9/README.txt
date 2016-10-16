Part I

Question 1:
python3 count_gv.py /home/assignments/assignment9/SNV_indel.biallelic.vcf /home/assignments/assignment9/sv.reclassed.filtered.vcf 

Variants:	Count
SNVs		5488226
INDELs		646978
DEL			1998
DUP			565
INV			106
MEI			1764
BND			2763.5
-
Question 2:
Proportion of genomic variants that are SVs: 0.001620786686897209
-
Question 3:
The spectrum of SVs in the individual NA12878 listed from greatest to least: BNDs, DELs, MEIs, DUPs, and INVs. I did not realize so much genetic variations could occur, especially the wide range of variant classes. The BND mutation count should be even given that for every BND, two sequences are reported. However, I got an odd count leading me to believe that at some point during pre-filtering, before I got the file, some of the BND mates were removed leaving only 1 behind. 
-
Question 4:
The vast majority of indels are 1 base in length and the number of occurances quickly decrease as the length increases, looks similar to a poisson distribution. I only plotted lengths 50 or less based on the definition of indels.

The vast majority of DELs are between 50 and 1,000 and the number of occurances quickly decreases as the length increases, looke similar to a poisson distribution as well. I only plotted lengths between 50 to 15,000 given that one or two DELs lengths were around 3x10^7 disrupting the view of the graph.

The vast majority of MEIs are around 500 in length approximately with a small peak occuring around 6,000. I only plotted lengths between less than 7,000 to make the graph more readable. 

I speculate that the length distribution might be much smaller if we limited the data to exonic indels only, large indels occurrig within the exonic region would most likely be lethal or cause major damage. So I would imagine indels would be smaller within the exonic region and much larger within the intronic regions where they would cause less trouble.
-
Part II

Question 5:
python3 quantify_genotype.py /home/assignments/assignment9/SNV_indel.biallelic.vcf 

Homozygous reference count:	3470043
Heterozygous count:	2857948
Homozygous alternate count: 1638628
Missing alleles count: 91342
-
Question 6:
The difference in the number of homozygous alternate (or non-reference homozygous) and heterozygous SNVs and indels does make biological sense because the population is extremely similar genetically therefore heterozygous alleles would be more than homozygous alternate because a mutational event within a single parent would have to occur to increase the number of possible heterozygous event inherited while at least a mutational event in both parents would have to occur for a homozygous alternate to be inherited.
-
Question 7:
*Restricted my analysis to the autosomes*
Number of variants that clearly violate the rules of Mendelian segregation: 145028
-
Question 8:
Four potential reasons that could explain the Mendlian violations
1) A de novo mutation could have occured within the child only
2) A sequencing error could have occured
3) A variant calling error could have occured
4) Variants with low quality could have made it through the filtering increasing the number of false positive variants
-
Question 9:
Number of variants that clearly violate the rules of Mendelian segregation: 26352
-
Comments:
I think everything went okay...
-
Suggestions:
...?
-

