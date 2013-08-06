###############################
## MyPDBPareser     - v 0.1  ##
## A. W. Covert III, Ph. D   ##
## All rights reserved       ##
###############################

from re import split
from os import path

#store and parse data out of a PDB file--currently only stores the raw dump and a list of residues
class cMyPDBParser(object):

    ### Function - cMyPDBParser::__init__
    ### Purpose  - Instasiate a solution object from the context file--ie the seed solution
    ### Input    - Filename of pdb file; we assume that this file is "clean"
    ### Output   - a cMyPDBParer Object
    def __init__ (self, fileName):
        if(path.exists(fileName)):
            #dictionary of residual values for translation purposes
            self.__resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
                             'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
                             'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
                             'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

            #self.__aaCodes = list(self.__resdict.viewvalues())
            self.__aaCodes = list(self.__resdict.values())

            self.__fileName = fileName

            #nab the file
            fp = open(self.__fileName,'r')
            self.__fileDump = fp.readlines()
            fp.close()
            
            #get the file length
            size = len(self.__fileDump)

            i = -1

            #split the atom vecotrs
            for i in range(0,size):
                if(self.__fileDump[i][0:4] == "ATOM"):
                    self.__fileDump[i] =  split('\s+',self.__fileDump[i])
                else:   #end of usefull stuff
                    break

            if(i < size):
                #delete the extra info in the PDB file
                del self.__fileDump[i:size]

            #print len(self.__fileDump)

            self.__residue = self.__parseResidue()
        else:
            print("ERROR: could not find file %s" % (fileName))
            quit()

    ### Function - cMyPDBParser::parseResidue
    ### Purpose  - Get a list of single letter residue codes
    ### Input    - a PDB file dump
    ### Output   - A list of single letter residue codes
    def __parseResidue(self):


        try:
            #grab the last residue position
            lastPos = self.__fileDump[-1][5]
        except IndexError:
            print self.__fileDump
            print self.__fileName

        #print len(self.__fileDump)
        
        currRes = int(self.__fileDump[0][5])

        #store the chainIDs one per chain present
        currChain = self.__fileDump[0][4]

        result = []
        result.append([])
        result[-1].append(self.__resdict[self.__fileDump[0][3]])

        keys = []
        keys.append(currChain)
        
        #Loop through each position
        #  Snag the residue
        #  Skip ahead the next residue
        for x in range(1, len(self.__fileDump)):
            if(currRes != int(self.__fileDump[x][5])):
                
                #update the index of the current residue
                currRes = int(self.__fileDump[x][5])

                if(self.__fileDump[x][4] != currChain):
                    currChain = self.__fileDump[x][4]
                    keys.append(currChain)
                    result.append([])
                    
                #add the residue to the list
                result[-1].append(self.__resdict[self.__fileDump[x][3]])

        return dict([(keys[x], result[x]) for x in range(0,len(result))])

    ### Function - cMyPDBParser::getResidue
    ### Purpose  - residue accessor function
    ### Input    - none
    ### Output   - Residue list
    def getResidue(self):
        return self.__residue

    ### Function - cMyPDBParser::getResidue
    ### Purpose  - residue accessor function
    ### Input    - none
    ### Output   - Residue list
    def getAAcodes(self):
        return self.__aaCodes

    
        
        
