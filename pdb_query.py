#!/usr/bin python
### HISTORY ###################################################################
# Version       Date    Coder   Comments
# 1.0   09.09.2016      khoulahan     Initial development
#
# This script queries PDB and outputs PBD ids that match query 
# The current query looks for the following: 
# Organism: E. Coli
# No ligand in structure
# Only protein in structure
# Molecular Weight between 25000-70000 Da
# Wildtype 
# Applies 90% similarity filter to ensure proteins with 90% similarity are 
# represented by only one structure
#
# Usage: python pdb_query.pdb > pbd_ids.txt
#
#
### PREAMBLE ######################################################################################
import urllib2
from optparse import OptionParser

### COMMAND LINE ARGUMENTS ########################################################################
parser = OptionParser()
parser.add_option("-q", "--query", dest="query", help="xml query", metavar="FILE")

(options, args) = parser.parse_args()

### SET QUERY #####################################################################################
url = 'http://www.rcsb.org/pdb/rest/search'

query_file = open(options.query, "r")
query_text = query_file.read()

### RUN QUERY AND RETURN RESULTS ##################################################################
req = urllib2.Request(url, data=query_text)

f = urllib2.urlopen(req)

for line in f:
    print line[0:4]
