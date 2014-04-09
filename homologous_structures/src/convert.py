import sys
fasta_aln = sys.argv[1]
phy_aln = sys.argv[2]

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

Convert(fasta_aln, phy_aln, 'fasta', 'phylip')

