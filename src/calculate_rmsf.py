#!/usr/bin/python

import sys, numpy

from prody import *

def loadFiles(pdb_name, dcd_name):
  ensemble = parseDCD(dcd_name)
  structure = parsePDB(pdb_name)

  ensemble.setCoords(structure)
  ensemble.setAtoms(structure.calpha)
  ensemble.iterpose()

  return (ensemble, structure)

def calculateRMSFs(ensemble, structure):

  rmsf_all_sites = []

  for i in range(0, len(structure.calpha)):
    print(i)
    ensemble.setAtoms(structure.select('ca resnum ' + str(i)))

    rmsf = ensemble.getRMSFs()
    rmsf_all_sites.append(numpy.mean(rmsf))

  return rmsf_all_sites

def generateRMSFFile(pdb_file, dcd_file, out_file):
  (ensemble, structure) = loadFiles(pdb_file, dcd_file)

  all_rmsfs = calculateRMSFs(ensemble, structure)

  output_handle = open(out_file, 'w')

  output_handle.write('rmsf\n')

  for i in all_rmsfs:
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

    generateRMSFFile( pdb_file, dcd_file, out_file )

if __name__ == "__main__":
  main()

