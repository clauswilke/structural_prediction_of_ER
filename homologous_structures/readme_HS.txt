<<<<<<< HEAD
Written by Dariya Sydykova, 4/10/14.

RMSF calculations for pdb homologous structures. 

Contents of directories:
	
rmsf/
	Contains the RMSF calculations for pdb structures with 5 or more sequences at 5% difference of aligned amino acids. RMSFs were calculated by src/get_rmsf.r using sequences in alignments/. The naming convention '*_HS.rmsf' here means RMSF values from Homologous Structures, mapped to the representative structure. Here * indicates the name of the representative structure.
	First column of data in these files 'res_num' is the residue number in the representative structure.
	Second column of data in these files 'res_name' is the residue name in the representative structure.
	Third column of data in these files 'rmsf' is the site-specific rmsf value measured using homologous 

alignments/
	Contains alignments of homologous pdb structures. Each folder in this directory contains:
	- *_all_relevan_seq.fa contains aligned sequences from blast_hits/cutoff_blast_hits/*_cutoff_blast_hits.csv
	- *_all_unique_seq.fa contains all unique sequences from *_all_relevan_seq.fa
	- *_10%_unique_seq.fa contains sequences with 10% difference in aligned amino acids from *_all_relevan_seq.fa
	- *_5%_unique_seq.fa contains sequences with 5% difference in aligned amino acids from *_all_relevan_seq.fa
	- *_2%_unique_seq.fa contains sequences with 2% difference in aligned amino acids from *_all_relevan_seq.fa
	Here * represents abbreviation of the viral protein name (ex: Hepatitis C Protease - HCP).

blast_hits/
	Contains BLAST output for all viral sequences in the table format. all_blast_hits contains raw tables that are reported by BLAST when a pdb sequence is blasted against the Protein Data Bank. cutoff_blast_hits contains altered all_blast_hits tables  with cutoffs imposed based on the identity and alignment length values. Detailed description of the cutoff imposed is contained in cutoff_blast_hits/cutoffs.txt manually for each pdb structure. 
	
pdb_availablity/
	Contains text files that track the number of structures present for BLAST outputs and assorted alignments. This folder also contains all_pdb_availability.txt that captures available structures for all the viral proteins analyzed. 
=======
This folder contains all the work related to RMSF derived from homologous crystal structures. This work was done by by Dariya Sydykova.

The subfolders contain the following information:
- alignments
    Sequence alignments of homologous structures at various
    sequence-divergence cutoffs.
    
- pdb_availability
    Summary tables from blast search.

- rmsf
    RMSF values calculated from structural alignments of homologous structures.

- trees
    Phylogenetic trees for the sequence alignments.

- blast_hits
    Raw blast results.
    
- src
    All analysis scripts used to carry out the work.
    
- weights
    Phylogenetic weights calculated from the trees in "trees."
>>>>>>> df15ba604d4799558a739f4f4a8a7cc7d430ff5a
