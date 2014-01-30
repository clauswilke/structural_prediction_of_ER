alignments:

The folder contains alignments of pdb structures with 10%, 5%, and 2% differences in aligned amino acids. The folders are named "ten_unique_seq.fa", "five_unique_seq.fa", and "two_unique_seq.fa" for each protein respectively. The alignments also contain an alignment with all unique sequences, named "all_unique_seq.fa", which has sequences with difference in 1 amino acid. 

pdb_availability.csv:
Table that contain the number of all pdb structures, related pdb structures, all unique pdb structures, unique pdb structures with 105, 5%, and 2% differences for each protein. The $notes sections explains the way cutoff was performed. 

align_all_blasthits.r:
Takes in pdb structures that were picked from all BLAST hits (cutoff_BLAST_hits) and aligns them to each other. It also finds a reference sequence (refseq) and reference sequence id (refpdb) to use a first sequence that is used in comparison to everything else. 

sort_aln.r:
Sorts through the alignment build in align_all_blasthits.r to establish sequences with certain percentage difference. 
