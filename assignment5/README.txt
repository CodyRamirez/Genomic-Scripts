Part 1.0
python3 analyze_WGBS_methylation.py BGM_WGBS.bed
-
Question 1:
DNA methylation across chromosome 21 looks particularly high around 0.8 to 1 and has a left tail.
-
Question 2:
CpG coverage across chromosome 21 is highest to the left and tails off to the right.
-
Question 2.1:
0.08414904690309218 fraction of the CpGs have 0X coverage.
-
Part 1.1
bedtools intersect -a BGM_WGBS_CpG_methylation.bed -b CGI.bed -wa -wb | bedtools groupby -g 5,6,7 -c 4 -o mean > WGBS_CGI_methylation2.bed
-
Part 1.2
python3 analyze_CGI_methylation.py WGBS_CGI_methylation.bed
-
Question 3:
DNA methylation for CpGs in CGIs increase between 0.8 and 1.0 with average CGI methylation levels having a spike at 0. It almost looks identical to all the CpGs on chromosome 21's distribution of CpG methylation. 
-
Part 1.3.0
Gene promoters
python3 generate_promoters.py refGene.bed

Promoters are known to typically exist upstream of the gene and average roughly between 100-1,000 bp. Therefore, since promoter locations and sizes can vary I simple defined my promoters to exist 250 bp upstream of the gene start site and have a 500 bp length.
-
Promoter-CGI and non-promoter-CGI
Generating promoter-CGI:	bedtools intersect -a CGI.bed -b refGene_promoters.bed > promoters_CGI.bed
Generating non-promoter-CGI:bedtools intersect -v -a CGI.bed -b refGene_promoters.bed > non_promoters_CGI.bed

I merely used bedtools intersect, so if any regions overlapped between CGIs and refGene_promoters then I made them promoter-CGIs.
Non-promoter CGIs were any regions that did not overlap between CGIs and refGene_promoters.

Promoter-CGI average CpG methylation command: bedtools intersect -a BGM_WGBS_CpG_methylation.bed -b promoter_CGI.bed -wa -wb | bedtools groupby -g 5,6,7 -c 4 -o mean > average_promoter_CGI_methylation.bed
Non-promoter-CGI average methylation command: bedtools intersect -a BGM_WGBS_CpG_methylation.bed -b non_promoters_CGI.bed -wa -wb | bedtools groupby -g 5,6,7 -c 4 -o mean > average_non_promoter_CGI_methylation.bed

python3 analyze_CGI_methylation.py average_promoter_CGI_methylation.bed
python3 analyze_CGI_methylation.py average_non_promoter_CGI_methylation.bed
-
Question 4:
The DNA methylation profiles of promoter-CGIs and non-promoter-CGIs are completely opposite of each other. CpG Island promoters on average are not methylated, while CpG Island non-promoters on average are methylated.
-
Part 1.3.1
bedtools getfasta -fi hg19_chr21.fa -bed promoter_CGI.bed -fo promoter_CGI.fa
bedtools getfasta -fi hg19_chr21.fa -bed non_promoters_CGI.bed -fo non_promoters_CGI.fa

python3 nuc_count_multisequence_fasta.py promoter_CGI.fa
Dinucleotide Frequencies:	CG:0.10283709394470994

python3 nuc_count_multisequence_fasta.py non_promoters_CGI.fa
Dinucleotide Frequencies:	CG:0.0960872366990213
-
Question 5:
DNA methylation located in a gene promoter can typically repress gene transcription. However, DNA methylation is essential for normal development. 60 to 80 percent of all CpGs are methylated within human somatic cells. However, GC-rich regions that possess a high density of CpGs remain methylation free.
-
Part 2
perl bed_reads_RPKM.pl CGI.bed BGM_MeDIP.bed > MeDIP_CGI_RPKM.bed
perl bed_reads_RPKM.pl CGI.bed BGM_MRE.bed > MRE_CGI_RPKM.bed

python3 compare_methylome_technologies.py MeDIP_CGI_RPKM.bed MRE_CGI_RPKM.bed WGBS_CGI_methylation.bed

The correlation between MeDIP-seq RPKM and MRE-seq RPKM is:(Spearman’s (rho) rank correlation = -0.68475909098635979, p-value = 3.9363119725264965e-51)
The correlation between MeDIP-seq RPKM and WGBS methylation is:(Spearman’s (rho) rank correlation = 0.8014612059891536, p-value = 6.1848072534316125e-82)
The correlation between MRE-seq RPKM and WGBS methylation is:(Spearman’s (rho) rank correlation = -0.84054942176344005, p-value = 2.4848690869671547e-97)

I picked Spearman's rank-order correlation because measures the strength of association between two ranked variables. My variables are ratios and I assumed there would be a monotonic relationship between my variables.
-
Question 6:
MeDIP-seq and methylation have a positive correlation.
MRE-seq and methylation have a negative correlattion
MeDIP-seq and MRE-seq have a negative correlation.
-
Outliers
Location: chr21 	9825442		9826296
I did not look further into issue.

The correlation between MeDIP-seq RPKM and MRE-seq RPKM is:(Spearman’s (rho) rank correlation = -0.69921908891505424, p-value = 5.5798207106731297e-54)
The correlation between MeDIP-seq RPKM and WGBS methylation is:(Spearman’s (rho) rank correlation = 0.80452640616407733, p-value = 8.6929820394681309e-83)
The correlation between MRE-seq RPKM and WGBS methylation is:(Spearman’s (rho) rank correlation = -0.84356905689255612, p-value = 2.0206563488510149e-98)
-
Comments:
I think things went okay.
-
Suggestions:
Y'all are doing a great job.
-
Extra credit
{Commands for running bed_reads_RKPM.pl}
{Command for running analyze_H3K4me3_scores.py}
{Copy analyze_H3K4me3_scores.py, H3K4me3_RPKM_promoter_CGI.bed, H3K4me3_RPKM_non_promoter_CGI.bed, and H3K4me3_RPKM_promoter_CGI_and_H3K4me3_RPKM_non_promoter_CGI.png}
-
Question EC.1:
{How does the H3K4me3 signal differ in promoter-CGIs and non-promoter-CGIs?}
-
Question EC.2:
{What are some better alternatives to model MeDIP-seq data and MRE-seq data instead of using RPKM? Explain.}
-
Question EC.3:
{What would be a better way to compare H3K4me3 values instead of using boxplots? Explain.}
-
