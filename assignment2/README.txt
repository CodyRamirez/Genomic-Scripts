Part 1:
Question 1:
I used the non-redundant database to save time searching. This way, I only compare my sequence to all non-redundant sequences instead of aligning to multiply sequences that are redundant.
-
Question 2:
125
-
Question 3:
Naumovozyma castellii, 820, 52%
-
Question 4:
124; I got fewer hits than before because a higher BLOSUM index requires "closer" matches because the matrix is based on more closely related sequences, therefore narrowing the number of hits.
-
Question 5:
Torulaspora delbrueckii
-
Question 6:
881, 69%
-
Question 7:
118; I think the number of hits decreased because lowering the gap existence penalty to 7 also increased the gap exstension to 2, thus making any gaps larger than 5 base pairs would be penalized greater.
-
Question 8:
For example, the score of the closest ortholog decreased from 820 to 732 because the same gap penality goes from 19 to 23 this is due to the gap cost changing from (extistence:11 extension:1) to (existence:7 extension:2)
-
Question 9:
I would expect the search to take more time if I lowered the word length because I imagine it would take longer to align more smaller reads and more time to realize a given sequence has a lower score than another read, thus wasting time.
-
Question 10:
YES!
-


Part 2:
Question 11
bowtie2 -x chr22_index -U reads.fq -S chr22.sam 2> sum_report.txt


7,113


10,131


9,418


-
Question 12
The dataset is enriched with nucleotides A and T along with dinucleotides AA, AG, CA, CC, CT, GA, GG, TG, and TT.

Nucleotide Fold Enrichment:
A:1.0636
C:0.9324
G:0.9436
T:1.0596

Dinucleotide Fold Enrichment:
AA: 1.272
AC: 0.814656
AG: 1.197712
AT: 0.970288
CA: 1.220816
CC: 1.0697824
CG: 0.260176
CT: 1.180448
GA: 1.004848
GC: 0.867232
GG: 1.088832
GT: 0.81496
TA: 0.757008
TC: 0.979088
TG: 1.229248
TT: 1.272688

The enrichment scores were calculated by dividing observed by expected number of frequencies. For example, the observed nucleotide frequency was divided by the expected nucleotide frequency of 0.25, while dinucleotides were divided by their expected dinucleotide frequency of 0.0625.

The dataset was enriched for A and T nucleotides along with AA, AG, CA, CC, CT, GA, GG, TG, and TT dinucleotides which I found that  methylation assays, like bisulfite sequencing, and even in vitro nucleosome-reconstitution assays are typically enriched for these same nucleotides and dinucleotides from what I could find.
-
Extra Credit 1:
Species: Canis familiaris
We used BLASTn instead of BLASTx because we are merely comparing a nucleotide query sequence against a nucleotide sequence database to identify where the sequence came from or is related too; while BLASTx was used mainly to translate the sequence into a peptide and then blast it against a protein database.
-
Comments:
It was unclear to me whether enrichment refered to the nuc_count.py program or whether I was suppose to calculate enrichment another way. Also, I found a lot of assays existed that met the enrichment results so I just listed some.
-
Suggestions:
I just think some clarification is required. Some of the questions seemed ambiguous and open ended for interpretation.
-
