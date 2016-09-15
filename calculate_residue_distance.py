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
structure = parser.get_structure('Test', "/home/keh/PDB/ZevTest/VAT_Hexamer_refined_210316-updatedATP_SHEET_HELIX_ANNOTATED.pdb")
#structure = parser.get_structure(options.code, options.pdb)

### GET FASTA SEQUENCE ############################################################################
# get sequence from fasta file used to calculate EC couplings
for fasta_record in SeqIO.parse("/home/keh/PDB/ZevTest/VAT_sequence.fasta", "fasta"):
#for fasta_record in SeqIO.parse(options.fasta, "fasta"):
    fasta_sequence = fasta_record.seq

### CONVERT THREE LETTER CODE TO ONE LETTER #######################################################
d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
def shorten(x):
    if len(x) % 3 != 0: 
        raise ValueError('Input length should be a multiple of three')
    y = ''
    for i in range(len(x)/3):
            y += d[x[3*i:3*i+3]]
    return y

### GET PDB SEQUENCE AND ALIGN WITH FASTA #########################################################
# get sequence from pdb file
all_sequences = {}
score_max = 0
# go through each chain, extract the amino acid sequence and find chain with the best
# alignment to the fasta file
# only uses the first model (can update this if necessary)
chain = model[0]
for chain in model:
    residues = []
    for residue in chain:
        # get residue name, check that it is a valid amino acid and convert to single letter 
        # amino acid code, het atoms do not get included in sequence
        residue_code = residue.get_resname()
        if residue_code in d:
            residues.append(shorten(residue_code))
        else:
            continue
    # check is any residues
    if len(residues) > 0:
        chain_sequence = ''.join(residues)
        # align the chain sequence with the fasta (or full sequence), use penalties for gaps
        alignment = pairwise2.align.globalms(fasta_sequence, chain_sequence, 2, 0, -.5, -.1)
        # check is alignment is better than previous alignments 
        # keeps first alignment if there is a tie
        if alignment[0][2] > score_max:
            chain_id = chain.get_id()
            score_max = alignment[0][2]
            best_alignment = alignment
        all_sequences[chain.get_id()] = ''.join(residues)

## RE-INDEX FASTA TO MATCH PDB ALIGNMENT ##########################################################
# re-index fasta residue indexing according to alignment
alignment_num = 0
fasta_num = 0
pdb_num = 0
indexes = pd.DataFrame({'Fasta':np.repeat('-', len(str(best_alignment[0][0]))),
                        'PDB': np.repeat('-', len(str(best_alignment[0][1])))})
for f,p in zip(str(best_alignment[0][0]),str(best_alignment[0][1])):
    if f != '-':
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

