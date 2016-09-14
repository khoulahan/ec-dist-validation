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
import pandas as pd
import numpy as np

### COMMAND LINE ARGUMENTS ########################################################################
parser = OptionParser()
parser.add_option("-p", "--pdb", dest="pdb", help="path to pdb file", metavar="FILE")
parser.add_option("-f", "--fasta", dest="fasta", help="path to fasta file used to calculate EC scores", metavar="FILE")
parser.add_option("-c", "--code", dest="code", help="pdb code of structure", default = "0000")
parser.add_option("-e", "--ecscores", dest="ecscores", help="path to EC scores file", metavar="FILE")
parser.add_option("-a", "--atom", dest="atom", help="atom to calculate distance from", default = "CA")

(options, args) = parser.parse_args()

### DATA ANALYSIS #################################################################################
parser = PDB.PDBParser()

# Parse the structure into a PDB.Structure object
#structure = parser.get_structure('Test', "/home/keh/PDB/test_subset/Test.pdb")
structure = parser.get_structure(options.code, options.pdb)

### GET FASTA SEQUENCE ############################################################################
# get sequence from fasta file used to calculate EC couplings
#for fasta_record in SeqIO.parse("/home/keh/PDB/test_subset/Test.fasta", "fasta"):
for fasta_record in SeqIO.parse(options.fasta, "fasta"):
    fasta_sequence = fasta_record.seq

### GET PDB SEQUENCE AND ALIGN WITH FASTA #########################################################
# get sequence from pdb file 
score_max = 0
record_num = 0
#for pdb_record in SeqIO.parse("/home/keh/PDB/test_subset/Test.pdb", "pdb-seqres"):
for pdb_record in SeqIO.parse(options.pdb, "pdb-seqres"):
    pdb_sequence = pdb_record.seq
    # align sequence from pdb to fasta file sequence
    alignment = pairwise2.align.localxx(pdb_sequence, fasta_sequence, one_alignment_only=True)
    # check is score of alignment is greater than the previous sequence 
    # this only applies when there are multiple chains reported in a pdb file
    if alignment[0][2] > score_max:
        score = alignment[0][2]
        chain_num = record_num
        score_max = score
        best_alignment = alignment
    record_num += 1

## RE-INDEX FASTA TO MATCH PDB ALIGNMENT ##########################################################
# re-index fasta residue indexing according to alignment
alignment_num = 0
fasta_num = 0
pdb_num = 0
indexes = pd.DataFrame({'PDB': np.repeat('-', len(str(best_alignment[0][1]))),
                        'Fasta':np.repeat('-', len(str(best_alignment[0][0])))})
for p,f in zip(str(best_alignment[0][0]),str(best_alignment[0][1])):
    if c != '-':
        indexes.iat[alignment_num,0] = fasta_num 
        fasta_num += 1   
    if p != '-':
        indexes.iat[alignment_num,1] = pdb_num
        pdb_num += 1    
    alignment_num += 1

### GET EC SCORES #################################################################################        
# read in EC scores file 
#ecscores = open("/home/keh/PDB/test_subset/Test_CouplingScoresModified.csv", "r")
ecscores = open(options.ecscores, "r")

# open file for writing 
results_file = open(''.join([options.code, "_Distances.txt"]), "w")
missing_file = open(''.join([options.code, "_MissingResidues.txt"]), "w")

### EXTRACT CHAIN #################################################################################
# get chain corresponding to best alignment 
chains = structure.get_chains()
chain = chains.next()
for i in range(chain_num):
    chain = chains.next()

### CALCULATE DISTANCES ###########################################################################
# iterate over coupled residues
for line in ecscores:
    # find residue coordinates from ev coupling
    # minus one ev coupling index because python is 0 based
    residueIndex = line.split(",")
    if indexes.iloc[int(residueIndex[0])-1,1] != '-' and indexes.iloc[int(residueIndex[1])-1,1] != '-':
        residueOneIndex = int(indexes.iloc[int(residueIndex[0])-1,0])
        residueTwoIndex = int(indexes.iloc[int(residueIndex[1])-1,0])
    else:
        # if either residue is not present in pdb, write line to missing file
        missingFile.write(line+"\n")
        continue
    # extract EC score to output along with distance
    ECscore = residueIndex[2]
    # find the residues in the PDB file
    residueOnePDB = chain.child_list[residueOneIndex]
    residueTwoPDB = chain.child_list[residueTwoIndex]
    print residueOnePDB
    print residueTwoPDB
    # calculate alpha carbon distance
    betaDist = residueOnePDB[options.atom] - residueTwoPDB[options.atom]
    # write to file
    results = [str(residueOnePDB.id[1]), residueOnePDB.resname, str(residueTwoPDB.id[1]), residueTwoPDB.resname, str(ECscore), str(betaDist)]
    resultsLine = "\t".join(results)
    resultsFile.write(resultsLine+"\n")

# close file
missing_file.close()
results_file.close()

