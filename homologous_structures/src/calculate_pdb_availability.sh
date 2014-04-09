#!/bin/bash
## $1 = input pdb name

PDB_ID=`echo $1 | grep -o "^[a-zA-Z0-9]\{4\}"`
PDB_CHAIN=$2
PDB_VIRUS_NAME=$3

echo "PDB_ID: $PDB_ID"
echo "PDB_CHAIN: $PDB_CHAIN"
echo "PDB_VIRUS_NAME: $PDB_VIRUS_NAME"

##make a directory for pdb files and split by chain pdb files
if [ ! -d "../raw_pdbs" ]; then
    mkdir ../raw_pdbs
    mkdir ../raw_pdbs/split_chain
fi

if [ ! -d "../blast_hits" ]; then
    mkdir ../blast_hits/
    mkdir ../blast_hits/all_blast_hits
    mkdir ../blast_hits/cutoff_blast_hits
fi

##find all blasthits for a given pdb
Rscript find_relevant_seqs.r $PDB_ID $PDB_CHAIN ${PDB_VIRUS_NAME}_all_blast_hits.csv
 
##Get all homologous pdb structures for a pdb from all blasthits and align them
if [ ! -d "../alignments" ]; then
    mkdir ../alignments
fi

if [ ! -d "../alignments/$PDB_VIRUS_NAME" ]; then
    mkdir ../alignments/$PDB_VIRUS_NAME
fi

##Filter out sequences from ${PDB_VIRUS_NAME}_all_blast_hits.csv. This part was done manually using guidelines of <alignmentlength> >= 90% and <identity> >= 35%. Use ../pdb_availability/cutoffs.txt for exact cuttoff guidelines for each pdb structure.

##Download all necessary pdb files from ${PDB_VIRUS_NAME}_cutoff_blast_hits.csv, align them, and output a fasta file in ../alignments/$PDB_VIRUS_NAME/$PDB_VIRUS_NAME_all_relevant_seq.fa.
Rscript align_all_relevant_seqs.r $PDB_VIRUS_NAME ${PDB_VIRUS_NAME}_cutoff_blast_hits.csv

##Find unique homologous pdb structures from all blasthits alignment
#all unique
Rscript find_unique_seqs.r $PDB_ID $PDB_CHAIN $PDB_VIRUS_NAME ${PDB_VIRUS_NAME}_all_relevant_seq.fa all

#unique at 10%
Rscript find_unique_seqs.r $PDB_ID $PDB_CHAIN $PDB_VIRUS_NAME ${PDB_VIRUS_NAME}_all_relevant_seq.fa 10

#unique at 5%
Rscript find_unique_seqs.r $PDB_ID $PDB_CHAIN $PDB_VIRUS_NAME ${PDB_VIRUS_NAME}_all_relevant_seq.fa 5

#unique at 2%
Rscript find_unique_seqs.r $PDB_ID $PDB_CHAIN $PDB_VIRUS_NAME ${PDB_VIRUS_NAME}_all_relevant_seq.fa 2
