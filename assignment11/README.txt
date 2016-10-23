Part 1
perl /home/assignments/assignment11/SNAP.pl /home/assignments/assignment11/YLR211C.alignment1.fasta part1/
perl /home/assignments/assignment11/SNAP.pl /home/assignments/assignment11/YLR211C.alignment2.fasta part1/
-
Question 1:
I believe that alignment1 is based on the protein sequence and alignment2 is based on the nucleotide sequence. I noticed that alignment1's sequence is always in blocks of 3 bases even the gaps, but alignment2's sequence has alignments that are not in blocks of 3 bases hinting to me that alignment1 was the converted protein sequence. Also I noticed that the synonymous rate it much higher in alignment1 reinforcing that it is the protein sequence because if a protein alignment occurred then the alignment would account for proper protein alignment not taking into account which bases coded for that amino acid.

Alignment2 suggests this gene has been evolving under positive selection given that its dn/ds ratio is greater than 1. Homologous genes with a dn/ds ratio above 1 are evolving under positive selection, meaning that at least some of the mutations concerned must be advantageous.
-
Part 2
Command to run run_SNAP.py on all of the alignment files in the alignments subdirectory:

python3 run_SNAP.py /home/assignments/assignment11/alignments part2/ 2> alignments.err 
-
Question 2:
The average dn/ds ratio is: 0.0844027567021

3 genes have dn/ds > 1; These are the genes: 
YHR097C	1.1821
YER179W	1.0982
YML094W	3.4578

The best evolutionary explanation for why YML094W has a dn/ds > 1 is because the gene is undergoing positive selection, meaning that at least some of the mutations concerned must be advantageous. If the mutations are neutral or disadvantageous, the ratio would be in the range of 0 to 1.
-
Part 3
Command to run plot_gene_length_vs_dnds.py on the output files from Part 2:

python3 plot_gene_length_vs_dnds.py alignments.err alignments_all_dnds.txt
-
Question 3:
The relationship between a gene's dn/ds ratio and its length seems to be a negative correlation.

The most likely explanation for this relationship is due to the fact that as the gene length get's smaller any little change will have a larger affect size on the over evolution, if the mutation is advantageous or not. It would take more nonsynonymous mutations to occur within a larger gene to cause an affect, than it would take in a smaller gene.

I propose to take gene length into account and normalize the data to improve one's ability to detect positive selection in long genes.
-
Part 4
Command to run calc_average_go_dnds.py on the S. cerevisiae GFF and output file from Part 2:}

python3 calc_average_go_dnds.py /home/assignments/assignment11/saccharomyces_cerevisiae.gff alignments_all_dnds.txt 
-
Question 4:
GO ID:0034553/GO term: mitochondrial respiratory chain complex II assembly has the highest dn/ds ratio at 0.9213
YBR044C is annotated with this GO

GO ID:0004347/GO term: glucose-6-phosphate isomerase activity and GO ID:0051173/GO term: positive regulation of nitrogen compound metabolic process have the lowest dn/ds ratio at 0.0025
YBR196C and YDR075W are annotated with their respective GO terms

I arrived at my answer by using sort to find the highest and lowest average dn/ds ratio in average_go_dnds.txt with:
sort -k 2 average_go_dnds.txt | head -n 5
sort -k 2 average_go_dnds.txt | tail -n 5
-
Comments:
I had trouble parsing out and properly gathering the GO IDs and their corresponding genes in a efficient manor and strongly believe my code is a bit clunky...
-

