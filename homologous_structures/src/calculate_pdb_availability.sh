#!/bin/bash
## $1 = input pdb name

PDB_ID=`echo $1 | grep -o "^[a-zA-Z0-9]\{4\}"`
PDB_CHAIN=`echo $2 | grep -o "[ABCD]\{1\}"`
PDB_VIRUS_NAME=`echo $3`

echo "PDB_ID: $PDB_ID"
echo "PDB_CHAIN: $PDB_CHAIN"
echo "PDB_VIRUS_NAME: $PDB_VIRUS_NAME"

##Find all blasthits for a given pdb 
#Rscript find_relevant_seqs.r <pdb_id> <pdb_chaind> <all_blast_hits.csv to written>
#Rscript find_relevant_seqs.r $PDB_ID $PDB_CHAIN ${PDB_VIRUS_NAME}_all_blast_hits.csv
 
##Get all homologous pdb structures for a pdb from all blasthits and align them
#Rscript align_all_relevant_seqs.r <all_blast_hits.csv to be written> <reference_seq to be written>
#mkdir ~/Desktop/Research/structural_prediction_of_ER_homologous_structures/alignments/${PDB_VIRUS_NAME}
#Rscript align_all_relevant_seqs.r $PDB_VIRUS_NAME ${PDB_VIRUS_NAME}_cutoff_blast_hits.csv

##Find unique homologous pdb structures from all blasthits alignment
#all unique
Rscript find_unique_seqs.r $PDB_VIRUS_NAME ${PDB_VIRUS_NAME}_ref_seq.fa ${PDB_VIRUS_NAME}_all_relevant_seq.fa all

#unique at 10%
Rscript find_unique_seqs.r $PDB_VIRUS_NAME ${PDB_VIRUS_NAME}_ref_seq.fa ${PDB_VIRUS_NAME}_all_relevant_seq.fa 10

#unique at 5%
Rscript find_unique_seqs.r $PDB_VIRUS_NAME ${PDB_VIRUS_NAME}_ref_seq.fa ${PDB_VIRUS_NAME}_all_relevant_seq.fa 5

#unique at 2%
Rscript find_unique_seqs.r $PDB_VIRUS_NAME ${PDB_VIRUS_NAME}_ref_seq.fa ${PDB_VIRUS_NAME}_all_relevant_seq.fa 2
