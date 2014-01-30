import os, re, subprocess
from sys import argv 
from cMyPDBParser import *

#	Originally Written by Art Covert and modified/edited by Eleisha Jackson
#	Last edited by Eleisha Jackson on August 24, 2013
#	Description: This code renames the sequences from designed structures made in Rosetta
#	This is used when you used different runs to create the proteins and want to re-number them.

#directory containing output for all runs
sourceDir = "../"

sourcePrefix = "yeast_RSA-"

no_chain = False #used to be false
#PDB = "1RD8_X"
temp = 0.6

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
    
    print resultDir
    print searchDir
    
elif no_chain:
    resultDir = sourcePrefix+protein+"_"+chain+"-"+temp+"-*"
    searchDir = sourceDir + protein + "_" + chain
    print resultDir
    print searchDir

#name of the expected output file
outputFilename = "_0001_last_0001\.pdb"

#bins = 10

files = []

outfile = "sequences_"+protein+"_"+chain+"_"+temp+".csv"

#fp = open(outfile,"w")
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
#This numbers all the designed files starting from 0.
count = 0
for file in files:
	new_filename = protein  + "_" + chain + "_" + temp + "_" + str(count) + "_0001_last_0001.pdb" 
	nameString = "cp " + file + " " + "/Users/Eleisha/Documents/Wilke_Lab/Amir_ER_Project/1RD8_X/1RD8_X_renamed/" + new_filename	
	process = subprocess.Popen(nameString, shell = True, stdout = subprocess.PIPE)
	process.wait() 
	count = count + 1
	
