import re, sys, subprocess, os, fnmatch, shutil
from Bio import AlignIO

filedict=['', 'crimean_congo.Nucleoprotein', 'H1N1_HA', 'H1N1_NP', 'HCV_NS5B', 'marburg.vp35', 'riftvalley', 'westnile_chainb', 'hivpol', 'JEV', 'westNile_ns2b', 'dengue_ns3']

file=filedict[int(sys.argv[1])]
alnfile = 'aln_nuc_'+file+'.fasta'
treefile='tree_'+file+'.tre'
final_file=file+'.out'

file1='/home/sjs3495/HyPhy/alignments/'+alnfile
command1='cp '+file1+' .'
run1=subprocess.call(command1, shell=True)

### FIRST MAKE THE TREE

# convert to phylip
aln = AlignIO.read(alnfile, 'fasta')
AlignIO.write(aln, 'temp.phy', 'phylip')
makeTree = '/share/apps/raxmlHPC-7.0.4/bin/raxmlHPC -m GTRGAMMA -s temp.phy -n out'
print "treeing"
subprocess.call(makeTree, shell=True)

## save the tree just in case!!!!
shutil.move('RAxML_result.out', '/home/sjs3495/HyPhy/'+treefile)

outtree = open('RAxML_result.out', 'r')
treestring=outtree.read()
outtree.close()



aln = open(alnfile, 'r')
alnlines = aln.read()
aln.close()

## write the input file
hyin = open('temp.fasta', 'w')
hyin.write(alnlines+'\n'+treestring)
hyin.close()	

callHyphy = '/home/sjs3495/bin/bin/HYPHYMP REL_5cat.bf > '+final_file
print callHyphy
subprocess.call(callHyphy, shell=True)
shutil.copy(final_file, '/home/sjs3495/HyPhy/RELresults/')
shutil.copy(hyin, '/home/sjs3495/HyPhy/RELresults/')
