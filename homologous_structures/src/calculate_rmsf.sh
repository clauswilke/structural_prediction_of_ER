#!/bin/bash
## $1 = input pdb names

PDB_ID=`echo $1 | grep -o "^[a-zA-Z0-9]\{4\}"`
PDB_CHAIN=$2
PDB_VIRUS_NAME=$3

echo "PDB_ID: $PDB_ID"
echo "PDB_CHAIN: $PDB_CHAIN"
echo "PDB_VIRUS_NAME: $PDB_VIRUS_NAME"

##Convert .fasta files to .phy (phylip) format to run raxml
#python convert.py ../alignments/$PDB_VIRUS_NAME/${PDB_VIRUS_NAME}_5%_unique_seq.fa ../alignments/$PDB_VIRUS_NAME/${PDB_VIRUS_NAME}_5%_unique_seq.phy

if [ ! -d "../tree" ]; then
    mkdir ../tree
fi

if [ ! -d "../tree/$PDB_VIRUS_NAME" ]; then
    mkdir ../tree/$PDB_VIRUS_NAME
fi

##Get the tree
#raxmlHPC-PTHREADS-SSE3 -T 2 -m PROTGAMMAWAG -s ${PDB_VIRUS_NAME}_5%_unique_seq.phy -n ${PDB_VIRUS_NAME} -# 20

#if [ -f "../*.${PDB_VIRUS_NAME}.RUN.*" ]; then
#    mv ../*.${PDB_VIRUS_NAME}.RUN.* ../trees/$PDB_VIRUS_NAME
#    mv ../*.${PDB_VIRUS_NAME} ../trees/$PDB_VIRUS_NAME
#fi

if [ ! -d "../weights" ]; then
    mkdir ../weights
fi

##Get the weights
#Rscript get_weights.r $PDB_VIRUS_NAME

if [ ! -d "../rmsf" ]; then
    mkdir ../rmsf
fi

##Calculate rmsf
Rscript get_rmsf.r $PDB_ID $PDB_CHAIN $PDB_VIRUS_NAME ${PDB_ID}_${PDB_CHAIN}_HS.rmsf

