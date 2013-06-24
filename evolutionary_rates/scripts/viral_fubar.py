# 6/22/13. Runs FUBAR. Goes with viral_fubar.qsub (array job for each file in filedict)

import os, sys, subprocess, shutil

rundir='FUBARmaterials/'

filedict=['', 'H1N1_HA', 'H1N1_NP', 'HCV_NS5B', 'HIV_REV_SJS', 'JEV', 'WestNile_ns2b', 'dengue_ns3']

file=filedict[int(sys.argv[1])]
alnfile='aln_nuc_'+prefix+'.fasta'
treefile='tree_'+prefix+'.tre'
final_file='fubar_'+prefix+'.csv'

file1='/home/sjs3495/ViralProject_SJS/alignments/'+alnfile
command1='cp '+file1+' '+rundir
run1=subprocess.call(command1, shell=True)
file2='/home/sjs3495/ViralProject_SJS/trees/'+treefile
command2='cp '+file2+' '+rundir
run2=subprocess.call(command2, shell=True)	

os.chdir(rundir)

shutil.move(alnfile, 'temp.fasta')
shutil.move(treefile, 'tree.tre')				
			
cline='/home/sjs3495/bin/bin/HYPHYMP autoFUBAR.bf CPU=2'
runit=subprocess.call(cline, shell=True)
shutil.move('tree.tre.fubar.csv', '../'+final_file)























