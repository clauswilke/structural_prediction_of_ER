Written by Stephanie Spielman, 6/22/13.

PRELIMINARY results for evolutionary rate analysis of viral sequences*

Contents of directories:

scripts/
	alntree.py	  
		->	Culls sequences, aligns, and builds trees	
		->	Culling entails...removes duplicates, anything with a stop codon (this step can be improved for increased yield), any sequence with >=5% ambiguous characters.
	
	viral_fubar.py and viral_fubar.qsub
		-> Runs fubar for all sequences. Requires the fubar materials which are currently in Stephanie's home directory, but I'll share them as needed.


unprocessed/
	Contains totally raw (the downloaded sequences, nothing done to them) nucleotide data files for all viruses*

processed/
	Contains culled (in alntree.py) nucleotide and protein fasta files. Unaligned.

alignments/
	Contains all nucleotide and protein alignments in fasta format. Created by "mafft --auto" using protein sequences, then those were back-translated.

nucaln/
	Contains only the nucleotide alignments (same as the nucleotide alignments in "alignments/")

trees/
	Contains all phylogenies. Created by "FastTree -nosupport -gamma -mlacc 2 -slownni -spr 4".

fubar_results/
	Contains all .csv fubar output files

*NOTRUN/
	Contains the raw sequence file for RiftValley. This virus was not run because all sequences were out of frame. Can come back to it later.

####################################################################################
Number of sequences per virus, alphabetically:
Crimean	68
HCV	217
H1N1_HA	1200
H1N1_NP	958
HIV_REV	1331
JEV	146
Marburg	42
WestNile	209
*RiftValley - all sequences were out of frame, so this hasn't been run yet***
