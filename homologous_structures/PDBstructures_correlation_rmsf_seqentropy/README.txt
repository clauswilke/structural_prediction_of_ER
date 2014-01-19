Written by Daria Sydkova

This folder contains the correlations between PDB structure RMSF values and sequence entropy obtained from ~/evolutionary_rates/alignments. Each protein folder contains the following:
-alignment of protein structures with 5% difference ("*_final_aln.fa")
-weighted RMSF of these protein structures ("*_wtdRMSF.csv")
-weights that were used for RMSF calculations ("*_weights.csv")
-table where RMSF and sequence entropy values are aligned against a given protein ("*evrate_rmsf_corr_vals.csv")
-protein RMSF map ("*_rmsf_map.txt")
-protein sequence entropy map ("*_seqent_map.txt")
-pdb structure for mapping 

align_rmsf_seqent.r:
1. takes in fubar calculate evolutionary rate (~/evolutionary_rate/fubar_results/fubar_*), both protein maps, and the pdb structure.
2. makes a table with RMSF and sequence entropy values aligned based on the pdb structure.

get_rmsf.r:
1. takes in weights of the protein structures and their final alignment
2. structurally aligns proteins
3. calculates RMSF across each site
4. returns rmsf vector

get_seqent.r:
1. takes in an alignment in fasta format
2. calculates entropy for each column in the alignment 
3. returns seq_ent vector

