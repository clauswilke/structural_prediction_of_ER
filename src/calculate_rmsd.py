#!/usr/bin/python

import sys

from prody import *
from numpy import *

def loadFiles(pdb_name, dcd_name):
  dcd = Trajectory(dcd_name)
  structure = parsePDB(pdb_name)

  dcd.setCoords(structure)

  dcd.link(structure)

  return (dcd, structure)

def calculateRMSDs(dcd, structure, skip):

  rmsd_all_sites = zeros(len(structure.calpha))

  for site in range(1, len(rmsd_all_sites) + 1):
    print(site)
    
    site_positions = structure.select('resnum ' + str(site))
    
    locations = zeros(3 * (len(dcd) - skip)).reshape(len(dcd) - skip, 3)
    
    dcd.reset()
    for i in range(skip, len(dcd)):
      locations[i - skip] = calcCenter(site_positions)

    center = mean(locations, axis = 0)
    rmsd_all_sites[site - 1] = sqrt(sum((locations - center)**2))
   
  return rmsd_all_sites

def generateRMSDFile(pdb_file, dcd_file, out_file, skip):
  (dcd, structure) = loadFiles(pdb_file, dcd_file)

  all_rmsds = calculateRMSDs(dcd, structure, skip)

  output_handle = open(out_file, 'w')

  output_handle.write('rmsd\n')

  for i in all_rmsds:
    output_handle.write(str(i) + '\n')

  output_handle.close()

  return 0

def main():
  if len( sys.argv ) != 5:
    print '''

    You screwed up choosing input files.

    '''
    print "     ", sys.argv[0], "<pdb file> <dcd file> <skip> <output file>"
    
  else:
    pdb_file = sys.argv[1]
    dcd_file = sys.argv[2]
    skip     = int(sys.argv[3])
    out_file = sys.argv[4]

    generateRMSDFile( pdb_file, dcd_file, out_file, skip )

if __name__ == "__main__":
  main()

