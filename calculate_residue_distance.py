#!/usr/bin python
### HISTORY #######################################################################################
# Version       Date    Coder   Comments
# 1.0   09.09.2016      khoulahan     Initial development
#
# Calculates distance between CA from two residues as specified in inputted
# EV coupling scores
# NOTE: assumes only one sequence in fasta file, if more than one, will only
# look at last sequence 
#
#
### PREAMBLE ######################################################################################
from Bio import PDB, SeqIO, pairwise2
from optparse import OptionParser

### COMMAND LINE ARGUMENTS ########################################################################
parser = OptionParser()
parser.add_option("-p", "--pdb", dest="pdb", help="path to pdb file", metavar="FILE")
parser.add_option("-f", "--fasta", dest="fasta", help="path to fasta file used to calculate EC scores", metavar="FILE")
parser.add_option("-c", "--code", dest="code", help="pdb code of structure", default = "0000")
parser.add_option("-e", "--ecscores", dest="ecscores", help="path to EC scores file", metavar="FILE")
parser.add_option("-o", "--output", dest="output", help="path to output file", metavar="FILE")
parser.add_option("-a", "--atom", dest="atom", help="atom to calculate distance from", default = "CA")

(options, args) = parser.parse_args()

### DATA ANALYSIS #################################################################################
parser = PDB.PDBParser()

# Parse the structure into a PDB.Structure object
structure = parser.get_structure(options.code, options.pdb)

# get sequence from fasta file used to calculate EC couplings
for fasta_record in SeqIO.parse(options.fasta, "fasta"):
    fasta_sequence = fasta_record.seq

# get sequence from pdb file 
score_max = 0
record_num = 0
for pdb_record in SeqIO.parse(options.pdb, "pdb-seqres"):
    pdb_sequence = pdb_record.seq
    # check that pdb sequence equals fasta sequence
    if pdb_sequence == fasta_sequence:
        chain_num = 0
        offset = 0
        break
    else:
        # if not equal, align sequence from pdb to fasta file sequence
        alignment = pairwise2.align.localxx(pdb_sequence, fasta_sequence, one_alignment_only=True)
        # check is score of alignment is greater than the previous sequence 
        # this only applies when there are multiple chains reported in a pdb file
        if alignment[0][2] > score_max:
            score = alignment[0][2]
            offset = alignment[0][3]
            chain_num = record_num
            score_max = score
    record_num += 1
        
# read in EC scores file 
ecscores = open(options.ecscores, "r")

# open file for writing 
resultsFile = open(options.output, "w")

# get chain
chains = structure.get_chains()
for i in range(chain_num):
    chain = chains.next()

# iterate over coupled residues
for line in ecscores:
    # find residue coordinates from ev coupling
    # add offset to residue index in the case that the fasta sequence doesn't
    # perfectly match pdb sequence
    residueIndex = line.split(",")
    residueOneIndex = int(residueIndex[0])+offset
    residueTwoIndex = int(residueIndex[1])+offset
    # extract EC score to output along with distance
    ECscore = residueIndex[2]
    # find the residues in the PDB file
    residueOnePDB = chain.child_list[residueOneIndex]
    residueTwoPDB = chain.child_list[residueTwoIndex]
    # calculate alpha carbon distance
    betaDist = residueOnePDB[options.atom] - residueTwoPDB[options.atom]
    # write to file
    results = [str(residueOnePDB.id[1]), residueOnePDB.resname, str(residueTwoPDB.id[1]), residueTwoPDB.resname, str(ECscore), str(betaDist)]
    resultsLine = "\t".join(results)
    resultsFile.write(resultsLine+"\n")

# close file
resultsFile.close()

