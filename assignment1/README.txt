Please provide the exact command line arguments you used to generate your results.
python3 nuc_count.py hs_ref_GRCh38.p2_chr22.fa
python3 make_seq.py 1000000 0.27 0.23 0.24 0.26
-
Question 1:

Raw Counts
A: 11198038
C: 9819368
G: 9943117
T: 11155239
N: 53597
-
Question 2:

Nucleotide Frequencies
A: 0.266
C: 0.233
G: 0.236
T: 0.265
-
Question 3:

Dinucleotide Frequencies
AA: 0.080
AC: 0.051
AG: 0.075
AT: 0.061
CA: 0.076
CC: 0.067
CG: 0.016
CT: 0.074
GA: 0.063
GC: 0.054
GG: 0.068
GT: 0.051
TA: 0.047
TC: 0.061
TG: 0.077
TT: 0.079
-
Dinucleotide Frequencies
AA : 0.081
AC : 0.057
AG : 0.061
AT : 0.070
CA : 0.057
CC : 0.058
CG : 0.058
CT : 0.058
GA : 0.061
GC : 0.057
GG : 0.061
GT : 0.061
TA : 0.070
TC : 0.057
TG : 0.061
TT : 0.071
-
Compare the two lists of frequencies. What are the differences? Can you provide a biological explanation for these differences?:
The frequencies are nearly identical; however, the CG content in chr22 is noticeably different from the randomly generated sequence with the same nucleotide frequencies. A biological explanation for the difference is a phenomenon called CG suppression where CGs are the least frequent dinucleotide, making up less than 1% of all dinucleotides in humans.
-
Comments:
I do not feel anything went wrong.
-
Suggestions:
Dictionaries and how to utilize their simple commands.
-
