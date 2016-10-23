Usage:
python3 neutral_rate.py PRE1.aln 
-
Question 1:
Fraction of wobble positions conserved: 79/174 = 0.45
-
Question 2
Cutoff for 10bp sequence conserved >= 8
Explanation: I used python interactive mode and ran the following code to view the probability distribution of calling a 10 bp sequence significantly more conserved than expected under the neutral hypothesis at a p-value <= 0.05. I added the end values as long as the total was below 0.05. Therefore, 8 base pairs or more must be conserved for a 10 bp sequence to be significantly more conserved than expected.

Code I typed into interactive mode:
from scipy.stats import binom
for i in range(0,11):
	print( str(i) + '\t' + str(binom.pmf(i, 10, 0.45)))
-
Question 3:
Number of regions: 4
S_cer_consreved.txt output:
Start   Stop    Scer Sequence
530     547     CCGGGGCATCTTTCGCA
545     589     CAGCGAAATCTTTACGGTGGCAAAAAATAAAGAAAA-GTGAATA
582     594     GTGAATATTGAA
603     622     CAA-GAAAAGGGAGCACCT
-
Question 4:
The number of binding sites in conserved regions compared to the number in the entire promoter region varies for conserved regions. The first conserved region contains 2 out of 7 total binding sites. The second conserved region contain 3 out of 7 toatl binding sites. The third conserved region contains 0 out of 7 total binding sites. While, the fourth conserved region contains 1 out of 7 total binding sites. RPN4 TF fits with the annotation that PRE1 is involved in yeast proteasome assembly because it stimulates expression of proteasome genes. RPN4 is non-essential; null mutant has lower abundance of proteasome subunits and shows decreased proteolytic of reporters. I arrived at this answer using the http://www.yeastgenome.org/ and searching the functions of the TF I found from JASPAR.
-
Comments:
I think everything went okay... We shall see though...
-