######## Convert between sequence formats
#requires modules:
from Bio import SeqIO
from Bio.Alphabet import *
def Convert(infile, outfile, informat, outformat):

	in_handle=open(infile, 'r')
	in_parsed=list(SeqIO.parse(in_handle, str(informat)))
	out_handle=open(outfile, 'w')
	umm=SeqIO.write(in_parsed, outfile, str(outformat))
	in_handle.close()
	out_handle.close()
	
	return 0

#Convert('fasta_alignments/DPH_5%_unique_seq.fa', 'phy_alignments/DPH.phy', 'fasta', 'phylip')
#Convert('fasta_alignments/WNP_5%_unique_seq.fa', 'phy_alignments/WNP.phy', 'fasta', 'phylip')
#Convert('fasta_alignments/JEH_5%_unique_seq.fa', 'phy_alignments/JEH.phy', 'fasta', 'phylip')
#Convert('fasta_alignments/HCP_5%_unique_seq.fa', 'phy_alignments/HCP.phy', 'fasta', 'phylip')
#Convert('fasta_alignments/RVFN_5%_unique_seq.fa', 'phy_alignments/RVFN.phy', 'fasta', 'phylip')
Convert('fasta_alignments/HP_5%_unique_seq.fa', 'phy_alignments/HP.phy', 'fasta', 'phylip')

