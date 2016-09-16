#!/usr/bin python
### HISTORY #######################################################################################
# Version       Date    Coder   Comments
# 1.0   09.12.2016      khoulahan     Initial development
#
# Validation pipeline: 
# 1) Query PDB with user specified queary 
# 2) Download corresponding PDB and Fasta files for each ID in query output
# 3) Calculate EC score for each fasta file
# 4) Calculate distance between atoms with highest EC scores
#
#
### PREAMBLE ######################################################################################
#from optparse import OptionParser
from subprocess import call

### COMMAND LINE ARGUMENTS ########################################################################
parser = OptionParser()
parser.add_option("-q", "--query", dest="query", help="PDB database xml query", metavar="FILE")
parser.add_option("-p", "--pdb", dest="pdb", help="PDB file of protein of interest", metavar="FILE")
parser.add_option("-c", "--code", dest="code", help="pdb code of structure", default = "0000")
parser.add_option("-t", "--threshold", dest="threshold", help="EC score threshold")

(options, args) = parser.parse_args()

### PDB QUERY #####################################################################################
# open file handle  
pdb_id_file = open("pdb_ids.txt", "w")
# run query and write IDs to file
call(["python", "pdb_query.py", "-q", options.query], stdout=pdb_id_file)

# read pdb ids file 
pdb_ids = open("pdb_ids.txt","r")

for line in pdb_ids:
    # remove new line character 
    line = line.rstrip('\n')
### DOWNLOAD DATA #################################################################################
    # download fasta file
    fasta_file = ''.join([line, ".fasta"])
    fasta_url = ''.join(["http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=",line])
    call(["curl", "-o", fasta_file, fasta_url])
    # download PDB file
    pdb_file = ''.join([line, ".pdb"])
    pdb_url = ''.join(["http://www.rcsb.org/pdb/files/", line, ".pdb"])
    call(["wget", pdb_url])
### CALCULATE EC SCORES ###########################################################################
    # run jackhmmer and calculate scores
### CALCULATE DISTANCES ###########################################################################
    # generate output file name
    distance_file = ''.join([line, "_Distances.txt"])
    # calculate distances between residues with top EC scores
    call(["python", "calculate_residue_distance.py", "-p", pdb_file, "-f", fasta_file, "-c", line,
        "-e", ec_file, "-o", distance_file])

### PROTEIN OF INTEREST ###########################################################################
# extract sequence from protein of interest pdb file
call(["python", "convert_pdb_to_fasta.py", "-p", options.pdb, "-c", options.code]) 

### CALCULATE EC SCORES POI #######################################################################
# run jackhmmer and calculate scores
### CALCULATE DISTANCES POI #######################################################################
# generate file names
fasta_file_poi = ''.join([options.code, '.fasta'])
distance_file_poi = ''.join([options.code, '_Distances.txt'])
# calculate distances 
call(["python", "calculate_residue_distance.py", "-p", options.pdb, "-f", fasta_file_poi, "-c", 
    options.code, "-e", ec_file, "-o", distance_file_poi])

### VISUALIZE DATA ################################################################################
call(["Rscript", "compare_ec_distances.R", "-p", "pdb_ids.txt", "-t", options.threshold, 
    "-n", options.code])

