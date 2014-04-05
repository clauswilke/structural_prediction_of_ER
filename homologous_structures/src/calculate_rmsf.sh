#!/bin/bash
## $1 = input pdb names

PDB_ID=`echo $1 | grep -o "^[a-zA-Z0-9]\{4\}"`
PDB_CHAIN=`echo $2 | grep -o "[ABCD]\{1\}"`
PDB_VIRUS_NAME=`echo $3`

echo "PDB_ID: $PDB_ID"
echo "PDB_CHAIN: $PDB_CHAIN"
echo "PDB_VIRUS_NAME: $PDB_VIRUS_NAME"

##Convert .fasta files to .phy (phylip) format
#python convert.py

##Get the tree
#raxmlHPC-PTHREADS-SSE3 -T 2 -m PROTGAMMAWAG -s ${PDB_VIRUS_NAME}.phy -n ${PDB_VIRUS_NAME} -# 20

##Get the weights
#Rscript get_weights.r $PDB_VIRUS_NAME

##Calculate rmsf
Rscript get_rmsf.r $PDB_VIRUS_NAME  ${PDB_VIRUS_NAME}_wtd_rmsf.csv