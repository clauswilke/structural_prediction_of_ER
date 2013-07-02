#!/usr/bin/python

import sys, numpy

from prody import *

def loadFiles(pdb_name, dcd_name):
  ensemble = parseDCD(dcd_name)
  structure = parsePDB(pdb_name)

  ensemble.setCoords(structure)

  return (ensemble, structure)

def calculateRMSDs(ensemble, structure):

  rmsd_all_sites = []

  for i in range(0, len(structure.calpha)):
    print(i)
    ensemble.setAtoms(structure.select('ca resnum ' + str(i)))
    ensemble.superpose()

    rmsd = ensemble.getRMSDs()
    rmsd_all_sites.append(numpy.mean(rmsd))

  return rmsf_all_sites

def generateRMSDFile(pdb_file, dcd_file, out_file):
  (ensemble, structure) = loadFiles(pdb_file, dcd_file)

  all_rmsds = calculateRMSDs(ensemble, structure)

  output_handle = open(out_file, 'w')

  output_handle.write('rmsd\n')

  for i in all_rmsds:
    output_handle.write(str(i) + '\n')

  output_handle.close()

  return 0

def main():
  if len( sys.argv ) != 4:
    print '''

    You screwed up choosing input files.

    '''
    print "     ", sys.argv[0], "<pdb file> <dcd file> <output file>"
    
  else:
    pdb_file = sys.argv[1]
    dcd_file = sys.argv[2]
    out_file = sys.argv[3]

    generateRMSDFile( pdb_file, dcd_file, out_file )

if __name__ == "__main__":
  main()

