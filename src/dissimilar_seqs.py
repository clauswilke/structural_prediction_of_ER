#!usr/local/bin/python
import string, re, gzip, urllib, shutil

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
import dissimilar_seqs

def readAlignment(filename):
    align = {}
    for seq_record in SeqIO.parse(filename,"fasta"):
        if seq_record.id in align:
            print "Warning: skipping duplicated entry", seq_record.id
        else:
            align[seq_record.id] = seq_record.seq
            #print len(seq_record.seq)
    return align
	
'''These two functions are both needed to get a dataset of dissimiar sequences from a starting dataset which is comprised of sequences. The main functions is get_dissimilar_sequences. You give it a list of sequences as the input. It returns all of the sequences that are at least 75 percent different from every other sequence. The ouput is a list of sequences from the orignal dataset that satisfy this condiition.'''


#This function takes a sequence and compares it to all the sequences in a list. It returns a boolean that is equal to True or False. If True then the sequnce is at least 75 percent similar to at least one sequence in the least.
def compare_sequences(seq, natural_seqs):
    cutoff = int(len(seq)*0.25) #The cutoff for sequence identity divergence
    sequence_similar = True
    for line in natural_seqs: #This compares the seq to all of the sequences in the list AA by AA
        dissimilar_count = 0
        for i in xrange(0, len(min([seq, line], key = len))): #Must take the length of the smaller seq being compared
            #print i
            natural_acid = seq[i] 
            other_acid = line[i]
            if (natural_acid != other_acid): #Checks if identity is the same at that site. 
                dissimilar_count = dissimilar_count + 1
                if(dissimilar_count>cutoff):
                    sequence_similar = False
                    return sequence_similar
                else:
                    continue
        if(dissimilar_count > cutoff):
            sequence_similar = False
    return sequence_similar 

#This is a function that takes a list of sequencs and a new list. This list is a list of all of the sequences in the orginal that are at most 25 sequence identity when compared to every other sequence in the list
'''def get_dissimilar_sequences(natural_seqs):
    dissimilar_natural_seqs = []
    sequence_is_similar = False
    for seq in natural_seqs:
        seq_copy = list(natural_seqs) #This copies a new list that can be compared
        seq_copy.remove(seq) #Must take out the sequence else it will compare to itself (will be always True)
        sequence_is_similar = compare_sequences(seq, seq_copy) #Determines whether it is similar
        if sequence_is_similar == False:
            dissimilar_natural_seqs.append(seq)
    return dissimilar_natural_seqs'''

def main():
    if len( sys.argv ) != 3:
        print '''
  This program maps an amino-acid alignment to a chain in a pdb file.
  The output is a map that relates each column in the input alignment, where
  possible, to a position in the PDB chain. The program requires the alignment
  software MAFFT to be installed to run properly.

  Usage:'''
        print "     ", sys.argv[0], "<input alignment> <output reduced alignment>"
    else:
        filename = sys.argv[1]
        #filename = "aln_aa_H1N1_HA.fasta"
        outfilename = sys.argv[2]
	outfile = open(outfilename,'w')
	align = readAlignment(filename)
	these_keys = []
	counter1 = 0
	counter2 = 0
	for i in align:
		counter1 += 1
		this_copy = align.copy()
		this_copy.pop(i)
		if (compare_sequences(align[i],this_copy.values())):
			counter2 += 1
			print i,counter1,counter2
			these_keys.append(i)
			outfile.write('>' + "%s\n" % i)
			outfile.write("%s\n" % align[i])
	print these_keys
	#print 'Amir'
	#for item in outfile:
	#	outfile.write("%s\n" % item)
        return 0

if __name__ == "__main__":
  main()
