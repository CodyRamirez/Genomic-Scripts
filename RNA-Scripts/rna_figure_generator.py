#!/usr/bin/env python3

"""
DESCRIPTION
	This script generates RNA-Seq metric figures using GMS generated data as input

	If you need to generate the input data required for this program refer to these GMS scripts:
		genome model rna-seq rna-seq-metrics --help
		genome model rna-seq splice-junction-metrics --help


USAGE
	python3 rna_figure_generator.py [--NormTransCov] <NormalizedTranscriptCoverage.tsv> [--RnaSeqMetrics] <RnaSeqMetrics.tsv> [--SpliceJunctionMetrics] <SpliceJunctionMetrics.tsv> [--Split] <Position>
	
	Note: Order matters, the options specified should always be followed by their corresponding required data.

	Example:
		python3 rna_figure_generator.py --NormTransCov NormalizedTranscriptCoverage.tsv --RnaSeqMetrics RnaSeqMetrics.tsv --SpliceJunctionMetrics SpliceJunctionMetrics.tsv

OPTIONAL INPUTS
	NormTransCov
		Requires: NormalizedTranscriptCoverage.tsv

	RnaSeqMetrics
		Requires: RnaSeqMetrics.tsv

	SpliceJunctionMetrics
		Requires: SpliceJunctionMetrics.tsv

	Split
		Requires: The position of where the data should be split. The position is inclusive.


NormTransCov OUTPUT
	Figure 1) rna_transcript_coverage.svg
		A line graph of normalized coverage across transcripts
		x-axis = Normalized Position Across Transcript (5' -> 3')
		y-axis = Normalized coverage
		legend = Samples

RnaSeqMetrics OUTPUT
	Figure 1) rna_genomic_distribution_of_aligned_reads.svg
		A stacked bar graph of the distribution of read alignments relative to genomic features, including coding, intergenic, intronic, ribosomal, and UTRs
		x-axis = Samples
		y-axis = Percentage
		legend = Genomic features

	Figure 2) rna_total_and_aligned_reads.svg
		A stacked bar chart of total reads generated and aligned from RNA-Seq
		x-axis = Samples
		y-axis = Total Reads (Millions)
		legend = type of reads

SpliceJunctionMetrics OUTPUT
	Figure 1) rna_depth_of_gene_coverage.svg
		A clustered bar chart of the dynamic range of gene coverage at varying depths
		x-axis = Samples
		y-axis = Number of genes
		legend = Sequencing read depths
"""





import sys, argparse, numpy, pandas, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt





if(len(sys.argv) == 1):
	sys.exit("ERROR: YOU DID NOT INCLUDE DATA FOR FIGURE GENERATION\n" + __doc__)

#Initializing parser to argparse.ArgumentParser
parser = argparse.ArgumentParser()

#Specifying optional arguments
parser.add_argument("--NormTransCov", help = "Requires full file path to NormalizedTranscriptCoverage.tsv", required = False)
parser.add_argument("--RnaSeqMetrics", help = "Requires full file path to RnaSeqMetrics.tsv", required = False)
parser.add_argument("--SpliceJunctionMetrics", help = "Requires full file path to SpliceJunctionMetrics.tsv", required = False)
parser.add_argument("--Split", help = "Requires the position of where the data should be split", required = False, type = int)
parser.add_argument("--DetailedHelpDoc", help = "You have to also type in the the magic word after you request the detailed documentation option", required = False)

args = parser.parse_args()





def generate_rna_transcript_coverage_figure(DataFrame, DataSetName):
	fig = plt.figure()
	DataFrame.plot.line(title = 'Transcript Coverage', ax = fig.gca())
	plt.legend(loc = 'upper center', bbox_to_anchor = (0.5, -0.2), ncol = 5, columnspacing = 1.0, labelspacing = 0.0, handletextpad = 0.0, handlelength = 1.5, fancybox = True)
	plt.xlabel('Normalized Position Across Transcript\n(5\' -> 3\')')
	plt.ylabel('Normalized Coverage')
	fig.savefig('rna_transcript_coverage_'+ DataSetName +'.png', bbox_inches = 'tight')

def generate_rna_genomic_distribution_of_aligned_reads_figure(DataFrame, DataSetName):
	fig = plt.figure()
	DataFrame.plot.bar(title = 'Genomic Distribution of Aligned Reads', ax = fig.gca(), stacked = True)
	plt.legend(labels = ['Coding', 'Intergenic', 'Intronic', 'Ribosomal', 'UTR'], loc = 'upper center', bbox_to_anchor = (0.5, -0.45), ncol = 3, columnspacing = 1.0, labelspacing = 0.0, handletextpad = 0.0, handlelength = 1.5, fancybox = True)
	plt.xlabel('')
	plt.ylabel('Percentage (%)')
	fig.savefig('rna_genomic_distribution_of_aligned_reads_' + DataSetName + '.png', bbox_inches = 'tight')

def generate_rna_total_and_aligned_reads_figure(DataFrame, DataSetName):
	fig = plt.figure()
	DataFrame.plot.bar(title = 'Total Reads Generated and Aligned', ax = fig.gca(), stacked = True, legend = 'reverse')
	plt.xlabel('')
	plt.ylabel('Total Reads (Millions)')
	fig.savefig('rna_total_and_aligned_reads_ ' + DataSetName + '.png', bbox_inches = 'tight')

def generate_rna_depth_of_gene_coverage_figure(DataFrame, DataSetName):
	fig = plt.figure()
	DataFrame.plot.bar(title = 'Depth of Gene Coverage', ax = fig.gca())
	plt.legend(labels = ['1X', '10X', '30X', '100X'], loc = 'upper center', bbox_to_anchor = (0.5, -0.45), ncol = 4, columnspacing = 1.0, labelspacing = 0.0, handletextpad = 0.0, handlelength = 1.5, fancybox = True)
	plt.xlabel('')
	plt.ylabel('Number of Genes')
	fig.savefig('rna_depth_of_gene_coverage_' + DataSetName + '.png', bbox_inches = 'tight')

def generate_rna_depth_of_splice_junction_coverage_figure(DataFrame, DataSetName):
	fig = plt.figure()
	DataFrame.plot.bar(title = 'Depth of Gene Coverage with >= 50% Junctions Detected', ax = fig.gca())
	plt.legend(labels = ['1X', '2X', '10X', '20X', '50X', '100X', '500X', '1000X'], loc = 'upper center', bbox_to_anchor = (0.5, -0.45), ncol = 4, columnspacing = 1.0, labelspacing = 0.0, handletextpad = 0.0, handlelength = 1.5, fancybox = True)
	plt.xlabel('')
	plt.ylabel('Number of Genes')
	fig.savefig('rna_depth_of_splice_junction_coverage_' + DataSetName + '.png', bbox_inches = 'tight')





if args.NormTransCov:

	NormTransCov_DataFrame = pandas.read_csv(args.NormTransCov, sep = '\t', index_col = 'POSITION')

	if args.Split:

		generate_rna_transcript_coverage_figure(NormTransCov_DataFrame.iloc[ : , :args.Split], "subset_1")
		generate_rna_transcript_coverage_figure(NormTransCov_DataFrame.iloc[ : , args.Split: ], "subset_2")

	else:

		generate_rna_transcript_coverage_figure(NormTransCov_DataFrame, "all_data")

	print("NormTransCov figure generation complete")





if args.RnaSeqMetrics:

	RnaSeqMetrics_DataFrame = pandas.read_csv(args.RnaSeqMetrics, sep = '\t', index_col = 'LABEL')

	RnaGenomicDistributionOfAlignedReads_DataFrame = RnaSeqMetrics_DataFrame.loc[ : , ['PCT CODING BASES', 'PCT INTERGENIC BASES', 'PCT INTRONIC BASES', 'PCT RIBOSOMAL BASES', 'PCT UTR BASES']]

	RnaTotalAndAlignedReads_DataFrame = RnaSeqMetrics_DataFrame[['TOTAL READS']].sub(RnaSeqMetrics_DataFrame['TOTAL READS MAPPED'], axis = 0).join(RnaSeqMetrics_DataFrame['TOTAL READS MAPPED'])
	RnaTotalAndAlignedReads_DataFrame = RnaTotalAndAlignedReads_DataFrame.loc[ : , ['TOTAL READS MAPPED', 'TOTAL READS']]

	if args.Split:

		generate_rna_genomic_distribution_of_aligned_reads_figure(RnaGenomicDistributionOfAlignedReads_DataFrame.iloc[ :args.Split, : ], "subset_1")
		generate_rna_genomic_distribution_of_aligned_reads_figure(RnaGenomicDistributionOfAlignedReads_DataFrame.iloc[ args.Split: , : ], "subset_2")

		generate_rna_total_and_aligned_reads_figure(RnaTotalAndAlignedReads_DataFrame.iloc[ :args.Split, : ], "subset_1")
		generate_rna_total_and_aligned_reads_figure(RnaTotalAndAlignedReads_DataFrame.iloc[ args.Split: , : ], "subset_2")

	else:

		generate_rna_genomic_distribution_of_aligned_reads_figure(RnaGenomicDistributionOfAlignedReads_DataFrame, "all_data")

		generate_rna_total_and_aligned_reads_figure(RnaTotalAndAlignedReads_DataFrame, "all_data")

	print("RnaSeqMetrics figure generation complete")





if args.SpliceJunctionMetrics:


	RnaDepthOfGeneCoverage_DataFrame = pandas.read_csv(args.SpliceJunctionMetrics, sep = '\t', index_col = 'LABEL')

	GeneCoverage_DataFrame = RnaDepthOfGeneCoverage_DataFrame.iloc[ : , [8, 7, 9, 6]]

	SpliceJunctionCoverage_DataFrame = RnaDepthOfGeneCoverage_DataFrame.iloc[ : , [21, 23, 20, 22, 25, 19, 24, 18]]

	if args.Split:

		generate_rna_depth_of_gene_coverage_figure(GeneCoverage_DataFrame.iloc[ :args.Split, : ], "subset_1")
		generate_rna_depth_of_gene_coverage_figure(GeneCoverage_DataFrame.iloc[ args.Split: , : ], "subset_2")

		generate_rna_depth_of_splice_junction_coverage_figure(SpliceJunctionCoverage_DataFrame.iloc[ :args.Split, : ], "subset_1")
		generate_rna_depth_of_splice_junction_coverage_figure(SpliceJunctionCoverage_DataFrame.iloc[ args.Split: , : ], "subset_2")

	else:

		generate_rna_depth_of_gene_coverage_figure(GeneCoverage_DataFrame, "all_data")
		
		generate_rna_depth_of_splice_junction_coverage_figure(SpliceJunctionCoverage_DataFrame, "all_data")

	print("SpliceJunctionMetrics figure generation complete")





if args.DetailedHelpDoc:
	print(__doc__)




