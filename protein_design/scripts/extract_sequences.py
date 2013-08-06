import os, re, subprocess
from sys import argv 
from cMyPDBParser import *

#	Originally Written by Art Covert and modified/edited by Eleisha Jackson
#	Last edited by Eleisha Jackson on August 8, 2013
#	Description: This code extracts the sequences from designed structures made in Rosetta

#directory containing output for all runs
sourceDir = "../"

sourcePrefix = "yeast_RSA-"

no_chain = False #used to be false

#a hack if we didn't specify the chain to examine -- needed for Austin's flue BS
if len(argv) == 4:
    protein = argv[1]
    
    chain = argv[2]

    temp = argv[3]
elif len(argv) == 3:
    
    no_chain = True #used to be true
    
    protein = argv[1]
    
    chain = "X"

    temp = argv[2]
else:
    print "Wrong number of arguments."
    quit()

#search_string_for the results directories
if not no_chain:
    resultDir = sourcePrefix+protein+"_"+chain+"-"+temp+"-*"
    searchDir = sourceDir + protein + "_" + chain
    
    #print resultDir
    #print searchDir
    
elif no_chain:
    resultDir = sourcePrefix+protein+"_"+chain+"-"+temp+"-*"
    searchDir = sourceDir + protein + "_" + chain
	#print resultDir
	#print searchDir

#name of the expected output file
outputFilename = "_0001_last_0001\.pdb"

#bins = 10

files = []

outfile = "sequences_"+protein+"_"+chain+"_"+temp+".csv"

fp = open(outfile,"w")
#find all of the PDB files output by backrub
#start by walking throught the directory of search results
for path, names, filename in os.walk(searchDir,False):

    #is this a result directory
    sPath = re.search(resultDir,path)
    if(sPath != None):

        #print path, names, filename
        #is the output file there, if it is tack it onto the list
        for file in filename:
            #print "searching '" + file + "' for string '" + outputFilename + "'" 
            sPath = re.search(outputFilename,file)
            if(sPath != None):
                files.append(path+"/"+file)

#print files
print "Found " + str(len(files)) + " files for processing"

residues = []
for file in files:
    pdb = cMyPDBParser(file)
    str = ''
    keys = pdb.getResidue().keys()
    for key in keys:
        for res in pdb.getResidue()[key]:
            str = str + res
    fp.write(str+"\n")
    residues.append(str)
fp.close()

