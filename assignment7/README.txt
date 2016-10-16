Question 1:
Score for highest affinity binding site for the TF scoring matrix = 48
Sequences corresponding to highest affinity binding site for the TF scoring matrix:
GCTGCGCACG
GCTGTGCACG

Score for highest affinity binding site for the polymerase scoring matrix = 62
Sequence(s) corresponding to highest affinity binding site for the polymerase scoring matrix:
GCACGGCACG

I added up all of the highest numbers from each column and determined the base from which row the highest number was in.
Rows: 1)A, 2)C, 3)G, 4)T
Using the transcription factor scoring matrix: 5+6+5+5+3+4+4+3+7+6 which corresponded to two sequences.
Using the polymerase scoring matrix: 3.5+6+5+7+6+5.5++7+7+8+7 which corresponded to a single sequence.
-
Question 2:
COMMAND:
python3 scan_sequence.py tf_score_matrix.txt promoter1.txt 40
OUTPUT:
orientation     sequecne        position        score
forward         GCTATGCACG      68              45.00
orientation     sequence        position        score
reverse         TGGTCTGTGA      36              43.00

COMMAND:
python3 scan_sequence.py tf_score_matrix.txt promoter2.txt 40
OUTPUT:
orientation     sequecne        position        score
forward         GCTGCGCACG      181             48.00
No threshold-exceeding hits found in the reverse direction!

COMMAND:
python3 scan_sequence.py polymerase_scoring_matrix.txt promoter1.txt 45
OUTPUT:
orientation     sequecne        position        score
forward         CCACGGCACG      102             59.00
No threshold-exceeding hits found in the reverse direction!

COMMAND:
python3 scan_sequence.py polymerase_scoring_matrix.txt promoter2.txt 45
OUTPUT:
orientation     sequecne        position        score
forward         GCACGGCACG      186             62.00
No threshold-exceeding hits found in the reverse direction!

{Name of the promoter that you'd expect to be repressed by the TF}
{Name of the promoter that you'd expect to be activated by the TF}
{Explanation}
-
Comments:
{Things that went wrong or you can not figure out}
-
Suggestions:
{What programming and/or genomics topics should the TAs cover in the next class that would have made this assignment go smoother?}
