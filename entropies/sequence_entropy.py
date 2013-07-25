#!/usr/bin/python
# This is a novice python code that reads in a pdb file with a given pdb id and the corresponding fasta sequence alignment file.
# Then uses DSSP to calculate RSA for the pdb file, also entropy for the alignment.
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
#import dissimilar_seqs
#import map_alignment_to_structure


def readAlignment(filename):
    align = {}
    for seq_record in SeqIO.parse(filename,"fasta"):
        if seq_record.id in align:
            print "Warning: skipping duplicated entry", seq_record.id
        else:
            align[seq_record.id] = seq_record.seq
            #print len(seq_record.seq)
    return align

def seqent(align):
    seqs = align.values()
    #print align.values()
    L = len(seqs[0])	# Length of each sequence
    #print L
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
		#print sum_freq
		entropy = 0.
		for n in counts[i].values():
			#print n, sum_freq
			entropy += n*(math.log(float(n)/sum_freq))/sum_freq
		entropy = -entropy
		if entropy < 0.:
			print entropy
			raw_input ('Are you joking? Entorpy is Negative?')
			stop
		entropy_list.append(entropy/sum_freq)
    #print entropy_list
    #print len(entropy_list)
    #print counts
    #print len(counts)
    return entropy_list

	
#ELEISHA'S CODE:

#Example of use
#pdb_id = '1B4T'
#chain_id = 'A'
#[AA_List, RSA] = get_values(pdb_id, chain_id)

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

#This is a dictionary that has the amino acid for the key and the max solvent accessibility for this amino acid
#THIS HAS BEEN UPDATED. I AM USING THE NEW THEORETICAL NUMBERS FROM THE 2013 AUSTIN, STEPHANIE, MATT, WILKE PAPER. 
residue_max_acc = {'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, \
                   'C': 158.0, 'Q': 223.0, 'E': 224.0, 'G': 104.0, \
                   'H': 209.0, 'I': 197.0, 'L': 201.0, 'K': 237.0, \
                   'M': 218.0, 'F': 239.0, 'P': 159.0, 'S': 151.0, \
                   'T': 172.0, 'W': 282.0, 'Y': 263.0, 'V': 174.0}
#This calculates the RSA values for a PDB using DSSP. It returns a list of the Amino Acids and a list of their RSA values. 
"""def get_values(pdb_id):
    searchPDB = pdb_id + ".pdb"  #This is the pdb file that is parsered by dssp
    pdbLocation = "/home/elj299/RSA_from_ceres/other/structs_all/" + searchPDB #This is the location of the pdb file that dssp will be parsing. 
    outputFile = pdb_id + "_DSSP.txt" 
    #processString = 'dssp' + ' -i ' + '"' + pdbLocation +'"'  + ' -o pdbOutput.txt '
    processString = 'dssp' + ' -i ' + searchPDB + ' -o ' + outputFile
    process = subprocess.Popen(processString, shell = True, stdout = subprocess.PIPE)
    process.wait() # Wait until dssp is done processing the file and calculating the Solvent Acessiblility  values
    #input = open("pdbOutput.txt" , 'r')
    #print os.path.exists(outputFile)
    #raw_input('Press Enter')
    input = open(outputFile, 'r')    
    fileContents = input.readlines()	
    string = fileContents[25]
    SAValue1 = string[13]
    SAList = [] #This is is the list which will store the SA values for each site
    AAList = [] #This is the list which will store the amino acid values for each site
    index = 0
    NoRSA = 0
    for line in fileContents:
        if index<25: #This skips the first few lines that do not have the SA value data
            index = index + 1
            continue
        else:  #Goes through each line with has the SA for each amino acid in order
            string = line #This stores the current line in the string "string"
            SAValue = string[35:39] #This stores the SA value for the current 
            AA = string[13] #This stores what the amino acid type is at the current  position
            number = int(SAValue)  #This turns the string fir the SA into an int type
            if AA !=( '!' or '*'): #This takes out the missing gaps that dssp might put in
                max_acc = residue_max_acc[AA] #This uses the dictionary to find the max SA for the amino acid at the current position (site)
                SAList.append(number/max_acc) #This divides the SA value for that position by the amx SA value position for tha amino acid. This normalizes the values and gives us the Relative Solvent Accessability (RSA) value. We this appends this value to the list of RSA values
                AAList.append(AA) #This appends the amino acid to the list
            else:
           	NoRSA = NoRSA + 1 #Counts the number of residues that did not have SA values
	    index = index + 1
    input.close() #Close the file with the dssp output file
    #os.remove('pdbOutput.txt') #Deletes the dssp output file 
    return (AAList, SAList) #Return the RSA values and the SAList"""

def main():
	#print sys.argv[0]
	#align_fname = sys.argv[1]
	#pdb_fname = sys.argv[2]
	#pdb_chain = sys.argv[3]
	#map_out_fname = sys.argv[4]
	#mapAlignmentToPDB(align_fname,pdb_fname,map_out_fname,pdb_chain)
	filename = "aln_aa_H1N1_HA.fasta"	# Note that this 
	align = readAlignment(filename)
	#print align
	#print align['CY039909']
	#reduced_align = dissimilar_seqs.get_dissimilar_sequences(align)
	#print reduced_align
	#print len(align), len(align)
	#raw_input('Amir')
	entropy = seqent(align)
	#print entropy
	pdb_id = '1RD8_I'
	#chain_id = 'A'
	#[AA_List, RSA] = get_values(pdb_id)
	#print len(RSA)
	#[AA_List, RSA] = get_values(pdb_id)
	outfilename = filename + '.entropy'
	outputfile=open(outfilename,'w')
	for item in entropy:
		outputfile.write("%s\n" % item)
	return

main()


# raw_input("Enter sequence filename: ")
#seqfile = open(filename,'r')
"""def readAlignment(filename):
	align = {}
	#sequences = []
	print os.path.exists(filename)
	for seq_record in SeqIO.parse(filename, 'fasta'):
		sequences.append(numpy.array(list(str(seq_record.seq))))
		
		#print seq_record.id
		if seq_record.id in align:
			print "Warning: skipping duplicated entry", seq_record.id
		else:
			align[seq_record.id] = [seq_record.seq]
			#print seq_record.seq
			
#	for position in range(0, len(sequences[1])):
#		print sequences[position]
		
	return sequences"""