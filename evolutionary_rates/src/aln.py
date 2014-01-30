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

def removeSeqs(raw_seqfile, nuc_raw, aa_raw, mapnames_file):
	'''Reads in an unaligned fasta file of nucleotide sequences. Removes all sequences which are shorter than 90% of the maximum sequence length as well as all sequences with any ambiguous nucleotides. Also removes all duplicate sequences.'''
	'''Returns unaligned nucleotide and protein sequence files.'''
	
	# First, removes sequences with too many ambiguities or that are <80% max sequence length
	# Second, checks reading frame. If not in frame but could be in frame, adjusts the sequence.
	# Third, removes any duplicate sequences. 

	#Read in the sequence file. Find the longest sequence length. Remove any sequences which are <0.8 that length. Additionally remove any sequences which contain "N"
	records=list(SeqIO.parse(open(raw_seqfile, 'rU'), "fasta"))
	newrecords=[]
	raw_numseqs=len(records)
	print raw_numseqs

	#print "Removing sequences with too many ambiguities or which are too short"
	cull=getMaxLength(records, '0.8')
	ambig_nucs=['Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'N']
	records=cullAmbigShort(records, cull, ambig_nucs)
	print "cull ambig short:", len(records)
	
	#print "Assuring that all sequences are in frame"
	records=checkReadingFrame(records)
	print "Reading frame:", len(records)
	
	#print "Checking for duplicates"
	records=killDuplicates(records)
	print "nuc duplicates:", len(records)
	
	'''
	## Check for duplicates with the TRANSLATED sequences. I KNOW THAT THE FOLLOWING IS INEFFICIENT AND SLIGHTLY GROSS, BUT YOU KNOW WHAT, IT WORKS. :)
	aa_records=[]
	for record in records:
		aaseq=(record.seq).translate()
		aa_record=SeqRecord(Seq(str(aaseq), generic_protein), id=record.id, description='')
		aa_records.append(aa_record)
	badIDs=killDuplicatesID(aa_records)
	
	final_records=[]
	for entry in records:
		if entry.id in badIDs:
			continue
		else:
			final_records.append(entry)
	print "aa duplicates:", len(final_records)
	'''

	# Write processed nucleotide and amino acid sequences to files
	outmap = open(mapnames_file, 'w')
	outfile_nuc=open(nuc_raw, 'w')
	outfile_aa=open(aa_raw, 'w')
	newname=0	
	for entry in records:
		aaseq=(entry.seq).translate()
		outfile_nuc.write('>'+str(newname)+'\n'+str(entry.seq)+'\n')
		outfile_aa.write('>'+str(newname)+'\n'+str(aaseq)+'\n')
		outmap.write(str(newname)+'\t'+str(entry.id)+'\n')
		newname+=1
	outfile_nuc.close()
	outfile_aa.close()

	return 0

################################################################################################

## FUNCTIONS THAT OPERATE IN removeSeqs

def checkReadingFrame(records):
	'''The sequences aren't in reading frame.'''
	seqlist=[]

	for entry in records:
		finalseq=''
		for i in range(3):
			newSeq=entry.seq[i:]
			trans=newSeq.translate()
			numStop=trans.count('*')
			if numStop==0:
				finalseq=newSeq
				break
			elif numStop==1:
				if trans[-1]=='*':
					finalseq=newSeq
					break
			elif numStop>1:
				continue	
		if finalseq!='':
			#print newSeq
			newrecord=SeqRecord(newSeq, id=entry.id, description='')
			seqlist.append(newrecord)
	return seqlist


def killDuplicates(records):
	'''Remove duplicate sequences'''
	count=0
	for ref_record in records:
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
		count+=1
	return records

def killDuplicatesID(records):
	'''Returns a list of which id's should be removed. DOES NOT ACTUALLY REMOVE THE SEQUENCES.'''
	count=0
	killIDs=[]
	
	for ref_record in records:
		w_count = count+1
		while w_count < len(records):
			test_record = records[w_count]
			if ref_record.id != test_record.id:
				# If sequences are the same but different ids, keep only one of the records
				if str(ref_record.seq) == str(test_record.seq):
					killIDs.append(test_record.id)
					removed_rec = records.pop(w_count)
					w_count-=1
			w_count += 1
		count+=1
	return killIDs
	
	
		

def getMaxLength(records, cutoff):
	'''Remove sequences which are too short'''
	seqlengths=[]
	for entry in records:
		seqlengths.append(len(entry.seq))
	cull=float(cutoff)*(numpy.amax(seqlengths))
	return cull
	
def cullAmbigShort(records, cull, ambig_nucs):
	ok=0
	seqlist=[]
	# Check that the length is adequate and N's are absent. If so, write that sequences to the file.
	for entry in records:
		nucseq=entry.seq
		seqlength=len(nucseq)
		if seqlength >= cull:
			numAmbig=0
			for ambig in ambig_nucs:
				has_ambig=False
				numAmbig = nucseq.count(ambig)
				if numAmbig > 0:
					has_ambig=True
					break
			# Continue if ambiguities
			if has_ambig==False:
				seqlist.append(entry)
	return seqlist
	
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

def buildTree(nuc_alnfile, name, treefile):
	'''Builds a tree with RAxML using nucleotide (ugh, fine) sequences, with no support values because I don't think HyPhy/FUBAR likes those.'''
	
	file = AlignIO.read(nuc_alnfile, 'fasta')
	AlignIO.write(file, 'temp.phy', 'phylip-relaxed')

	# single raxml inference
	runRaxml = 'raxmlHPC -m GTRGAMMA -s temp.phy -n out'
	runtree=subprocess.call(str(runRaxml), shell=True)
		
	shutil.move('RAxML_bestTree.out', treefile)
	
	# remove vomit	
	command = 'rm RAxML_*'
	subprocess.call(command, shell=True)
	
	return 0
################################################################################################
################################################################################################

#seqfiles=GetFiles('fasta', '../Thursday/unprocessed/')
#pdbseqs = list(SeqIO.parse(open('../Sunday/pdbseqs.fasta', 'rU'), 'fasta'))

#'
seqfiles=['riftvalley.fasta', 'hivpol.fasta', 'H1N1_HA.fasta', 'dengue_ns3.fasta', 'H1N1_NP.fasta', 'HCV_NS5B.fasta', 'marburg.vp35.fasta', 'JEV.fasta', 'westNile_ns2b.fasta', 'westnile_chainb.fasta', 'crimean_congo.Nucleoprotein.fasta']
for file in seqfiles:
	find = re.search('(.+)\.fasta', file)
	if find:
		name = find.group(1)
		full_raw = '../unprocessed/'+name+'.fasta'
		raw_aa = '../processed_8.7/'+'rawaa_'+name+'.fasta'
		raw_nuc = '../processed_8.7/'+'rawnuc_'+name+'.fasta'
		aa_alnfile = '../alignments_8.7/'+'aln_aa_'+name+'.fasta'
		nuc_alnfile = '../alignments_8.7/'+'aln_nuc_'+name+'.fasta'
	
		## since many names may be too long for raxml etc to deal with in phylip format, replace with ints but save the map in case need to know later
		mapnames_file = '../processed_8.7/'+name+'_taxonmap.txt'
		
		print file
		removeSeqs(full_raw, raw_nuc, raw_aa, mapnames_file)
						
		print "Aligning"
		alignSeqs(raw_aa, raw_nuc, aa_alnfile, nuc_alnfile)

