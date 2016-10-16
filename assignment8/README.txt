Question 1:
68
-
Question 2:
Explanation of each MetaGeneMark flag:
-a Shows the protein sequence of predicted genes
-d Shows the nucleotide sequence of predicted genes
-f G Outputs GFF2 format
-m [filename] File with gene finding parameters
-o [filename] Outputs file name

I would have not added or dropped any flags for my particular problem or use.
-
Question 3:
19
-
Question 4:
Number of ORFs shared between MetaGeneMark and call_orfs.py output: 8
Number of ORFs unique to call_orfs.py output: 60
Number of ORFs unique to MetaGeneMark output: 11
Only a small portion of the ORFs are shared due the parameters of needing ORFs to be >100bp in the call_orfs.py since so few ORFs were output from MetaGeneMark and even fewer had lengths >100bp.
-
Question 5:
{qseqid} : Means Query Seq-ID
{sseqid} : Means Subject Seq-ID
{pident} : Mean Percentage of identical matches
{length} : Means Alignment lenth
{mismatch} : Means Number of mismatches
{gapopen} : Means Number of gap openings
{qstart} : Means Start of alignment in query
{qend} : Means End of alignment in query
{sstart} : Means Start of alignment in subject
{send} : Means End of alignment in subject
{evalue} : Means Expect value
{bitscore} : Means Bit score
{slen} : Means Subject sequence length
{stitle} : Means Subject Title
-
Question 6:
Number of genes in BLAST output: 493 
Number of genes after filtering in the Python script: 3 
-
Question 7:
HMMER is can be used to search sequence databases for homologs of protein sequences, such as searching a protein sequence database with a single protein profile HMM.
First, you will need to build a profile HMM with hmmbuild then search the sequence database with hmmsearch.

HMMER is can be used to search sequence databases for homologs of DNA sequences, such as searching a DNA sequence database.
First, you you can optionally build a profile HMM with hmmbuild then search the DNA sequence database with nhmmer.
-
Question 8:
Total number of genes annotated by Resfams: 7
Resfams is extremely less than BLAST, however after filtering the BLAST genes Resfams and BLAST are much closer in comparison, only differing by a total of 4 genes.
-
Comments:
I had a problem figuring out how to append or grab all of the sequence from the fasta formatted file: mgm_orfs.faa
I just briefly made my mgm_orfs.faa sequence exist on a single line so I could run my program.
-
Suggestions:
{What programming and/or genomics topics should the TAs cover in the next class that would have made this assignment go smoother?}


