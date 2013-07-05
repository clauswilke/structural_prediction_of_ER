#!/usr/bin/python
import sys, subprocess, os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.Polypeptide import *
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import Dice

def parsePDBStructure( pdb_fname ):
    parser = PDBParser()
    structure = parser.get_structure('test', pdb_fname)
    return structure


def renumberResidues( structure, new_pdb ):
    model = structure[0]

    i = 1
    for chain in model:
        for residue in chain:
            residue.id = (' ', i, ' ')
            i += 1

    w = PDBIO()
    w.set_structure(structure)
    w.save(new_pdb)

    return structure


def getPDBSequence( pdb_fname, chain_id ):
    # first we create a renumbered structure of just the chain we're interested in
    structure = parsePDBStructure( pdb_fname )
    tmp_file = "tmp_sdielccskeoksmndivls.pdb"
    Dice.extract( structure, chain_id, 0, 100000, tmp_file )
    
    pdb_basename = os.path.basename( os.path.splitext( pdb_fname )[0] )
    new_pdb = pdb_basename + ".renumbered.pdb"
    structure = parsePDBStructure( tmp_file )
    renumbered_structure = renumberResidues( structure, new_pdb )
    
    polypeptides = ""
    ppb=PPBuilder()
    for pp in ppb.build_peptides(structure):
        polypeptides += pp.get_sequence()
    os.unlink( tmp_file )
    return polypeptides


def readAlignment( filename ):
    align = {}
    for seq_record in SeqIO.parse( filename, "fasta"):
        if seq_record.id in align:
            print "Warning: skipping duplicated entry", seq_record.id
        else:
            align[seq_record.id] = seq_record.seq
    return align


def createMafftInput( align_fname, pdb_record, mafft_in_fname ):
    # first write out pdb sequence
    out_handle = open( mafft_in_fname, "w" )
    SeqIO.write( [pdb_record], out_handle, "fasta" )

    # now just copy over the input alignment
    in_handle = open( align_fname, "r" )
    for line in in_handle:
        out_handle.write( line )

    out_handle.close()
    in_handle.close()

def matchSequences( seq_orig, seq_new, seq_pdb ):
    #print seq_orig
    #print seq_new
    #print seq_pdb
    
    #print "# columns in original alignment:", len( seq_orig )
    #print "# of amino acids in pdb reference:", len( str(seq_pdb).replace('-', '') )
    
    # first we test that the original and the new sequence are actually the same
    if str(seq_orig).replace('-', '') != str(seq_new).replace('-', ''):
        print "sequence in old and new alignment don't agree, can't match!"
        sys.exit()
    
    i = -1 # counts positions in old alignment
    j = -1 # counts positions in new alignment
    k = -1 # counts positions in pdb sequence
    match_list = []
    while True:
        i += 1
        # if we've reached the end of the original alignment we're done
        if i == len(seq_orig):
            break
        # if we're in a gap position at the old sequence, nothing needs to be done
        if seq_orig[i] == '-':
            #print i
            match_list.append( '-' )
            continue
        # we now have to increment j to match the same position, incrementing k
        # accordingly
        while True:
            j += 1
            # record whether this column in new alignment corresponds to an
            # occupied position in the pdb sequence
            if seq_pdb[j] != '-':
                k += 1
                pdb_match = True
            else:
                pdb_match = False
            if seq_new[j] != '-':
                assert seq_new[j] == seq_orig[i] # make sure numbering is still in sync
                break
        if pdb_match:
            #print "%i -> %i" % ( i, k )
            match_list.append( k )
        else:
            #print "%i -> NA" % i
            match_list.append( 'NA' )
    return match_list

def updateMatchCounts( pdb_match_counts, match_list ):
    '''Builds a histogram of which sites in the pdb are mapped how often
    to which sites in the original alignment'''
    i = 0
    for x in match_list:
        if x in pdb_match_counts:
            if i in pdb_match_counts[x]:
                pdb_match_counts[x][i]+=1
            else:
                pdb_match_counts[x][i]=1
        else:
            pdb_match_counts[x]={i:1}
        i += 1

def getConsensus( align ):
    '''Calculate the consensus sequence of an alignment. Gap characters can be the
    consensus if there are more gaps than any amino acid.'''
    seqs = align.values()
    L = len( seqs[0] )
    counts = [{} for i in xrange(L)]
    for seq in seqs:
        for i in xrange( L ):
            c = seq[i]
            if c in counts[i]:
                counts[i][c] += 1
            else:
                counts[i][c] = 1
    consensus=''
    for d in counts:
        sorted_counts = sorted( d.items(), key=lambda x: x[1], reverse=True )
        consensus += sorted_counts[0][0]
#    print consensus
    return consensus
    
        
        
def matchAlignments( align_fname, mafft_out_fname, pdb_seq_id, out_handle ):
    align_orig = readAlignment( align_fname )
    align_new = readAlignment( mafft_out_fname )
    pdb_match_counts = {}
    for key in align_orig:
        match_list = matchSequences( align_orig[key], align_new[key], align_new[pdb_seq_id] )
        updateMatchCounts( pdb_match_counts, match_list )
        
    # The variable pdb_match_counts now contains the full histogram of site counts.
    # All we have to do is extract the most common matches
    map_to_structure = {}
    for item in pdb_match_counts.items():
        pdb_position = item[0]
        sorted_counts = sorted( item[1].items(), key=lambda x: x[1], reverse=True )
        if pdb_position != '-' and pdb_position != 'NA':
            map_to_structure[sorted_counts[0][0]] = pdb_position
    
    # output results
    pdb_seq = str(align_new[pdb_seq_id]).replace('-', '')
#    out_handle.write( str(pdb_seq) )
#    out_handle.write( "\n" )
    out_handle.write( "align_pos\tconsensus_aa\tpdb_pos\tpdb_aa\n" )
    consensus = getConsensus( align_orig )
    for i in xrange( len( consensus ) ):
        if i in map_to_structure:
            pdb_position = map_to_structure[i]
            pdb_char = pdb_seq[pdb_position]
            pdb_position += 1 # change numbering from python to sane
        else:
            pdb_position = "NA"
            pdb_char = "-"
        out_handle.write( "%i\t%s\t%s\t%s\n" % ( i + 1, consensus[i], str(pdb_position), pdb_char ) )

def mapAlignmentToPDB( align_fname, pdb_fname, map_out_fname, pdb_chain ):

    mafft_in_fname = "tmp_92lesklwekd_mafft.in"
    mafft_out_fname = "tmp_92lesklwekd_mafft.out"


    pdb_seq_id = pdb_fname

    align_orig = readAlignment( align_fname )
    pdb_sequence = getPDBSequence( pdb_fname, pdb_chain )
    pdb_record = SeqRecord( pdb_sequence, id=pdb_seq_id, description="" )
    createMafftInput( align_fname, pdb_record, mafft_in_fname )

    print "*************************************************************"
    print "Aligning sequences with pdb structure\n"
    command = "mafft --auto " + mafft_in_fname + " > " + mafft_out_fname
    subprocess.call(command, shell=True)

    print "*************************************************************"
    print "Building map\n"
    handle = open( map_out_fname, "w" )
    matchAlignments( align_fname, mafft_out_fname, pdb_seq_id, handle )
    handle.close()
    
    print "*************************************************************"
    print "Map written to", map_out_fname
    
    # cleanup
    os.unlink( mafft_in_fname )
    os.unlink( mafft_out_fname )

def main():
    if len( sys.argv ) != 5:
        print '''
  This program maps an amino-acid alignment to a chain in a pdb file.
  The output is a map that relates each column in the input alignment, where
  possible, to a position in the PDB chain. The program requires the alignment
  software MAFFT to be installed to run properly.

  Usage:'''
        print "     ", sys.argv[0], "<input alignment> <pdb file> <pdb chain> <output map file>"
    else:
        align_fname = sys.argv[1]
        pdb_fname = sys.argv[2]
        pdb_chain = sys.argv[3]
        if pdb_chain == '-': # use hyphen to indicate missing chain id
		pdb_chain = ' '
        map_out_fname = sys.argv[4]
        mapAlignmentToPDB( align_fname, pdb_fname, map_out_fname, pdb_chain )

main()
