#!/usr/bin python
### HISTORY #######################################################################################
# Version       Date    Coder   Comments
# 1.0   09.09.2016      khoulahan     Initial development
#
# Convert pdb to fasta
#
#
### PREAMBLE ######################################################################################
from Bio import SeqIO
from optparse import OptionParser

### COMMAND LINE ARGUMENTS ########################################################################
parser = OptionParser()
parser.add_option("-p", "--pdb", dest="pdb", help="path to pdb file", metavar="FILE")
parser.add_option("-c", "--code", dest="code", help="pdb code of structure", default = "0000")

(options, args) = parser.parse_args()

### DATA ANALYSIS #################################################################################
# set up input handle
pdb_handle = open(options.pdb, "rU")
# set up fasta file name and handle
fasta_file = ''.join([options.code, ".fasta"])
fasta_handle = open(fasta_file, "w")
# parse pdb
sequence = SeqIO.parse(pdb_handle, "pdb-seqres")
# write as fasta
SeqIO.write(sequence, fasta_handle, "fasta")

fasta_handle.close()
pdb_handle.close()
