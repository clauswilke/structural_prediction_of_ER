Written by Stephanie Spielman, 8/9/13.

Results for evolutionary rate analysis of viral sequences

Contents of directories, in order of importance:

alignments/
	Contains all nucleotide and protein alignments in fasta format. Created by "mafft --auto" using protein sequences, then those were back-translated.

siterates/
	Contains the site omegas for each protein, as calculated by 5 category REL model in HyPhy. Each tab-delimited-file shows the protein position/residue and its associated omega.

trees/
	Contains all phylogenies. Created by RAxML. Single GTRGAMMA inference.
	ADDITIONALLY contains files for phylogenetic weights calculated by BranchManager.

HyPhy/
	hyin/
		Contains input files given to hyphy (fasta alignment and the tree)

	hyout/
		Contains raw hyphy output files. Parse these with parseHyPhy.py
	
	HyPhyFiles/
		Contains all the files necessary for running HyPhy (model file, a file with some important functions, and the batch file with 5 rate categories)
		

scripts/
	aln.py	  
		-> Culls sequences and aligns. Alignment done with amino acid sequences and then back-translated to codon alignments.
		-> Culling entails...removes nucleotide duplicates and sequences with any ambiguous characters. Additionally places everything into correct reading frame and removes anything that can't be.
	
	buildTrees.py
		-> Small script to build trees with RAxML. A single inference on the nucleotide data is performed with model GTRGAMMA, on cluster!

	runHyPhy.py
		-> Implements REL model with 5 rate categories (estimated from data)
	
	parseHyPhy.py
		-> Deals with the hyphy output files. Will output a text file of two columns: position  omega

unprocessed/
	Contains totally raw (the downloaded sequences, nothing done to them) nucleotide data files for all viruses

processed/
	Contains culled nucleotide and protein fasta files. Unaligned. Culling is performed by aln.py.
	Additionally contains "taxonmap" files. Basically, the taxon names within the alignments are quite long. Since they must be converted into phylip format for RAxML, names will be truncated which can be problematic as some could end up the same. So I changed all taxon names to ints and mapped out the actual taxon names here just in case!!
	
