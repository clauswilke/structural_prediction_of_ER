#!/usr/bin/python
# This is a novice python code that reads in a pdb file with a given pdb id and the corresponding fasta sequence alignment file.
import sys, subprocess, os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.Polypeptide import *
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import Dice
import math
import numpy


def readAlignment(filename):
    align = {}
    for seq_record in SeqIO.parse(filename,"fasta"):
        if seq_record.id in align:
            print "Warning: skipping duplicated entry", seq_record.id
        else:
            align[seq_record.id] = seq_record.seq
    return align

def seqent(align):
    seqs = align.values()
    L = len(seqs[0])	# Length of each sequence
    counts = [{} for i in xrange(L)]
    entropy_list = []
    for seq in seqs:
        for i in xrange(L):
            c = seq[i]
            if c not in "ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy":
				continue
            if c in counts[i]:
                counts[i][c] += 1
            else:
                counts[i][c] = 1
	entropy = 0.
    for i in xrange(L):
		sum_freq = sum(counts[i].values())
		entropy = 0.
		for n in counts[i].values():
			entropy += n*(math.log(float(n)/sum_freq))/sum_freq
		entropy = -entropy
		if entropy < 0.:
			print entropy
			raw_input ('Are you joking? Entorpy is Negative?')
			stop
		entropy_list.append(entropy/sum_freq)
    return entropy_list

def main():
	filename = "aln_aa_dengue_ns3.fasta"
	align = readAlignment(filename)
	entropy = seqent(align)
	outfilename = filename + '.entropy'
	outputfile=open(outfilename,'w')
	i = 0
	outputfile.write("res_num entropy\n")
	for item in entropy:
		i += 1
		outputfile.write(str(i) + " " + str(item) + "\n")
	return

main()