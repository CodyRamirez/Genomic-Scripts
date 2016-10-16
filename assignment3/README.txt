USAGE:
python3 map_sequence_starter.py /home/assignments/assignment3/scott_mouse_cDNAs.fa /home/assignments/assignment3/mckinley_raw_reads.txt
-
Question 1:
Read in the cDNAs
Created dictionary
Name	Reads	Reads per BP
Abcg2	7	0.0035460992907801418
Atp1A1	370	0.12044270833333333
Cd34	7	0.006092254134029591
ChAT	0	0.0
Gap43	36	0.05263157894736842
Gfap	38	0.029389017788089715
Mbp	304	0.4037184594953519
Myod1	0	0.0
Olig2	6	0.006172839506172839
Tubb3	154	0.11382113821138211

Question 1.2:
A.1)	Atp1A1

A.2)	The protein encoded by this gene belongs to the family of P-type cation transport ATPases. Na+/K+ -ATPase is an integral membrane protein responsible for establishing and maintaining the electrochemical gradients of Na and K ions across the plasma membrane. These gradients are essential for osmoregulation, for sodium-coupled transport of a variety of organic and inorganic molecules, and for electrical excitability of nerve and muscle.

A.3)	I think this gene is highly expressed in the sample because it is essential for electrical excitability of nerves which makes sense since the sample is mouse brain. However, the brain is made up of many cells, including neurons; neurons are cells that send and receive electro-chemical signals to and from the brain and nervous system.

B.1)	Mbp

B.2)	The protein encoded by Mbp gene is a major constituent of the myelin sheath of oligodendrocytes and Schwann cells in the nervous system.

B.3)	I think this gene is highly expressed in the sample because it is essential for myelin sheath of oligodendrocytes and Schwann cells in the nervous system which are major parts of the brain. 

C.1)	Tubb3

C.2)	This gene encodes a class III member of the beta tubulin protein family. Beta tubulins are one of two core protein families (alpha and beta tubulins) that heterodimerize and assemble to form microtubules. This protein is primarily expressed in neurons and may be involved in neurogenesis and axon guidance and maintenance.

C.3)	I think this gene is highly expressed in the sample because it is essential for axon guidance and maintenance which are necessary for the brain to function properly. Axons are long slender projections of a nueron that conducts electrical impulses away from the neuron's cell body. Axons are in effect the primary transmission lines of the nervous system and as bundles they make up nerves.

Question 1.3:
A.1)	Abcg2

A.2)	This gene encodes proteins that transport various molecules across extra- and intra-cellular membranes. This protein functions as a xenobiotic transporter which may play a major role in multi-drug resistance. Xenobiotic transport that may play an important role in the exculsion of xenobiotics from the brain.

A.3)	I think this gene is lowely expressed because it is an over all transporter protein that's main function as a xenobiotic transporter and if there are low levels of xenobiotics within the brain at the time of sample, then the expression level of the gene would also be low.

B.1)	ChAT

B.2)	This gene codes for a protein called choline acetyltransferase. This protein is located at the ends of nerve cells in specialized areas called presynaptic terminals. Choline acetyltransferase facilitates the production of a molecule called acetylcholine which is essential for normal muscle movement.

B.3)	I think this gene has basically no expression at all in the brain because it is essential for normal muscle movement and located at the ends of nerve cells in the presynaptic terminals so I would imagine even though they are receiving their signals from the brain, expression would not exist in brain tissue.

C.1)	Myod1

C.2)	This gene encodes for a nuclear protein that belongs to the basic helix-loop-helix family of transcription factors and the myogenic factors subfamily. It regulates muscle cell differentiation by inducing cell cycle arrest, a prerequisite for myogenic initiation. The protein is also involved in muscle regeneration. It activates its own transcription which may stabilize commitment to myogenesis.

C.3)	I think this gene has basically no expression at all in the sample because its a protein associated with muscle cells and the sample is from a brain. Therefore, I would not expect any expression of this gene within a brain tissue sample.

-
Question 2:
Yes; There is an enrichment of genes annotated as having "brain specific" expression in their list of the top 20 most expressed genes

Chi-Square Statistical Test; I used this test because there are two categorical variables from a single population and this test is used to determine whether there are significant association between the two variables. Therefore, I am determining whether the tissue sample is related to gene expression.

Two-tailed

	   No Expression | Expression |
Not Brain |	260	 |     10     | 270
Brain     |	20       |     10     |  30
	  	280      |     20     | 300

Expected Cell Frequencies:
	  |	252	|   18    |
	  |	28      |    2    |

(260-252)^2/252 + (20-28)^2/28 + (10-18)^2/18 + (10-2)^2/2 = Chi-square statistic 38.0952381

p-value = 0.00000206; p-value is < 0.00001; result significant at p < 0.05
-
Question 3:
We calculated hits as the raw number of counts divided by length of the gene, rather than just using raw number of counts because the relative expression of a transcript is proportional to the number of cDNA reads that originate from it. This helps to normalize the data. So if two genes have the same number of reads aligned to each of them, but one gene is half the length of the other then the smaller gene is more highly expressed than the larger gene. So the size of the gene and the number of reads have to be used to properly calculate expression or reads per base pair.
-
Question 4:
I am not sure what I would do exactly, but I know that I would first address the problem of redundant 25mer reads to ensure my dictionary did not become populated with repeatative reads effecting my search time. 

Limitation 1: Hash tables are fast in searching BUT it requires a perfect match, therefore making it insufficient when there are many collisions hurting alignment accuracy.  
Limitation 2: Hash tables take constant time on average BUT the cost of a good hash function can be significantly higher than the inner loop of the lookup algorithm for a sequential search.
Limitation 3: Hash tables can require lots of memory to store the indexed reads, thus slowing everything down slightly by having to search through the entire dictionary. Also repeatative 25mers could be created within the 25mer_dictionary causing time to be lost to search over non-redundant 25mers which could share multiply homologies with several genes, which then causes an assignment problem for which gene would receive a read count.
-
Comments:
I am really fuzzy on the limitations and how to address better programming methods for this assignment, basically question 4, I answered as best I could but I do not truly understand as much as I would like.
-
Suggestions:
I am not sure.

