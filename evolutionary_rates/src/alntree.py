import re, os, sys, subprocess, shutil, fnmatch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna,generic_protein
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
import numpy

################################################################################################
################################################################################################
def Pal2Nal(palfile, nucfile, paltype, nuctype, outfile, outputformat):

	aln_parsed=list(SeqIO.parse(str(palfile), str(paltype)))
	nuc_parsed=list(SeqIO.parse(open(nucfile, 'rU'), str(nuctype)))
	
	if len(aln_parsed)!=len(nuc_parsed):
		print palfile+' '+nucfile+' have different number of sequences! Please make sure that these two files correspond and that all stop codons are removed.'
		assert 1==0
	else:
		numseq=len(aln_parsed)
	
	nucMSA=MultipleSeqAlignment([])
	for p in range(0,numseq):
		pal_seq=str(aln_parsed[p].seq) #aa alignment sequence
		pal_id=str(aln_parsed[p].id)
		for n in range(0,numseq):
			if nuc_parsed[n].id==pal_id:
				nuc_seq=str(nuc_parsed[n].seq)
		nal=str()
		start=0 #counter for codon starting position
		end=3 #counter for codon ending position
		for position in pal_seq:
			#If gapped, missing, or masked position in alignment, append 3 gaps/missing/NNN to new string 
			if position=='-':
				nal=nal+'---'
			elif position=='?':
				nal=nal+'???'
			elif position=='X':
				nal=nal+'NNN'
			#If amino acid there, append corresponding codon
			else:
				codon=str(nuc_seq[start:end])
				nal=nal+codon
				start+=3
				end+=3
		#Make nucleotide MSA object
		aln_record=SeqRecord(Seq(str(nal), generic_dna), id=pal_id, description='')
		nucMSA.append(aln_record)
	#write alignment to file
	outfile=open(outfile, 'w')
	umm=AlignIO.write(nucMSA, outfile, str(outputformat))
	outfile.close()
	return 0
	

################################################################################################

def GetFiles(ext, dirpath):
	files=[]
	dir=os.listdir(dirpath)
	for file in dir:
		if fnmatch.fnmatch(file, '*.'+str(ext)):
			files.append(file)
	return files

################################################################################################

def removeSeqs(raw_seqfile, nuc_raw, aa_raw):
	'''Reads in an unaligned fasta file of nucleotide sequences. Removes all sequences which are shorter than 75% of the maximum sequence length as well as all sequences with any ambiguous nucleotides. Also removes all duplicate sequences.'''
	'''Returns unaligned nucleotide and protein sequence files.'''
	
	#Read in the sequence file. Find the longest sequence length. Remove any sequences which are <0.75 that length. Additionally remove any sequences which contain "N"
	records=list(SeqIO.parse(open(raw_seqfile, 'rU'), "fasta"))
	newrecords=[]
	raw_numseqs=len(records)
	
	## First, this removes all duplicate sequences and modifies records to contain only unique sequences and id's.
	print "Checking for duplicates"
	count=0
	for ref_record in records:
		removed = 0
		w_count = count+1
		while w_count < len(records):
			test_record = records[w_count]
			if ref_record.id != test_record.id:
				# If sequences are the same but different ids, keep only one of the records
				if str(ref_record.seq) == str(test_record.seq):
					removed_rec = records.pop(w_count)
					w_count-=1
			#If same ID and same seq, remove
			else:
				removed_rec = records.pop(w_count)
				w_count-=1
			w_count += 1
		count += 1
	# Get the acceptable length of sequence.
	seqlengths=[]
	for entry in records:
		seqlengths.append(len(entry.seq))
	cull=0.75*(numpy.amax(seqlengths))
	
	outfile_nuc=open(nuc_raw, 'w')
	outfile_aa=open(aa_raw, 'w')	
	
	print "Checking for lengths and N's"
	ok=0
	ambig_nucs=['Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'N']
	# Check that the length is adequate and N's are absent. If so, write that sequences to the file.
	for entry in records:
		nucseq=entry.seq
		seqlength=len(nucseq)
		maxAmbig=float(0.05*seqlength)
		if seqlength > cull:
			
			## Check for amibiguities
			numAmbig=0
			for ambig in ambig_nucs:
				numAmbig = nucseq.count(ambig)
			# Continue if ambiguities
			if numAmbig > maxAmbig:
				print "Too many ambiguities:", numAmbig
				continue
			
			# If last position is a stop codon, remove it
			lastcodon=str(nucseq[-3:])
			if lastcodon=='TAG' or lastcodon=='TAA' or lastcodon=='TGA':
				new_nucseq=nucseq[:-3]
			else:
				new_nucseq=nucseq
			
			# Find if any stop codons floating around and remove those sequences.
			aaseq=new_nucseq.translate()
			numStop = aaseq.count('*')
			if numStop==0:
				ok+=1
				outfile_nuc.write('>'+str(entry.id)+'\n'+str(new_nucseq)+'\n')
				outfile_aa.write('>'+str(entry.id)+'\n'+str(aaseq)+'\n')
			else:
				print "There were",numStop,"stop codons"
	outfile_nuc.close()
	outfile_aa.close()
	
	return (raw_numseqs, ok)

################################################################################################

def checkReadingFrame(records):
	'''The sequences aren't in reading frame.'''

	for entry in records:
		count=0
		seq = entry.seq
		print seq
		aaseq=seq.translate()
		print aaseq
		if aaseq[-1]=='*':
			aaseq=aaseq[:-1]
		numStop = aaseq.count('*')
		while numStop > 0:
			count+=1
			aaseq=aaseq[1:]
			if aaseq[-1]=='*':
				aaseq=aaseq[:-1]
			numStop = aaseq.count('*')
			if count > 2:
				print "big fail"
				break
		if count > 2:
			break
		else:
			print "achieved at count=",count

		return 0

################################################################################################

def alignSeqs(aa_new_seqfile, nuc_new_seqfile, aa_alnfile, nuc_alnfile):
	'''Generate an alignment from amino acid data. Return both an amino acid and nucleotide alignment (need Pal2Nal).'''
	
	# Make alignment
	alignCommand='mafft --auto --quiet '+aa_new_seqfile+' > '+aa_alnfile
	runalign=subprocess.call(alignCommand, shell=True)
	
	# Get the nucleotide alignment as well
	Pal2Nal(aa_alnfile, nuc_new_seqfile, 'fasta', 'fasta', nuc_alnfile, 'fasta')
	
	return 0
################################################################################################

def buildTree(aa_alnfile, treefile):
	'''Builds a thorough tree with FastTree using protein sequences, with no support values because I don't think HyPhy/FUBAR likes those.'''
	
	#Build Tree
	treeCommand='FastTree -nosupport -gamma -mlacc 2 -slownni -spr 4 '+aa_alnfile+' > '+treefile
	runtree=subprocess.call(treeCommand, shell=True)
	
	## Sometimes FastTree doesn't always make the tree, so make sure with this.
	size=os.path.getsize(treefile)
	while size==0:
		runtree=subprocess.call(treeCommand, shell=True)
		size=os.path.getsize(treefile)

	return 0
################################################################################################
################################################################################################

seqfiles=GetFiles('fasta', '../unprocessed/')

for file in seqfiles:
	print file
	find_name=re.search('(.+)\.fasta', file)
	if find_name:
		name = find_name.group(1)
		full_raw = '../unprocessed/'+file
		raw_aa = '../processed/'+'rawaa_'+name+'.fasta'
		raw_nuc = '../processed/'+'rawnuc_'+name+'.fasta'
		aa_alnfile = '../alignments/'+'aln_aa_'+name+'.fasta'
		nuc_alnfile = '../alignments/'+'aln_nuc_'+name+'.fasta'
		treefile = '../trees/'+'tree_'+name+'.tre'
		
		print "Culling sequences"
		(raw_numseqs, kept)=removeSeqs(full_raw, raw_nuc, raw_aa)
		
		print "Out of",raw_numseqs,", kept",kept
		print "Aligning"
		alignSeqs(raw_aa, raw_nuc, aa_alnfile, nuc_alnfile)
		print "Treeing"
		buildTree(aa_alnfile, treefile)

