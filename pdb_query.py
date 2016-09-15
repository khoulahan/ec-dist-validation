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
import urllib2

### SET QUERY #####################################################################################
url = 'http://www.rcsb.org/pdb/rest/search'

queryText = """

<orgPdbCompositeQuery version="1.0">
 <queryRefinement>
  <queryRefinementLevel>0</queryRefinementLevel>
  <orgPdbQuery>
    <queryType>org.pdb.query.simple.TreeEntityQuery</queryType>
    <description>TaxonomyTree Search for Escherichia coli (E. coli)</description>
    <t>1</t>
    <n>562</n>
    <nodeDesc>Escherichia coli (E. coli)</nodeDesc>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>1</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <queryType>org.pdb.query.simple.NoLigandQuery</queryType>
    <description>Ligand Search : Has free ligands=no</description>
    <haveLigands>no</haveLigands>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>2</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <version>head</version>
    <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
    <description>Chain Type: there is a Protein chain but not any DNA or RNA or Hybrid</description>
    <containsProtein>Y</containsProtein>
    <containsDna>N</containsDna>
    <containsRna>N</containsRna>
    <containsHybrid>N</containsHybrid>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>3</queryRefinementLevel>
  <orgPdbQuery>
   <queryType>org.pdb.query.simple.WildTypeProteinQuery</queryType>
   <description>Simple query for a list of PDB Entity IDs</description>
   <includeExprTag>Y</includeExprTag>
   <percentSeqAlignment>95</percentSeqAlignment>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>4</queryRefinementLevel>
  <orgPdbQuery>
   <queryType>org.pdb.query.simple.MolecularWeightQuery</queryType>
   <description>Molecular Weight Search : Min Molecular Weight=25000.0 Max Molecular Weight=70000.0</description>
   <mvStructure.structureMolecularWeight.min>25000.0</mvStructure.structureMolecularWeight.min>
   <mvStructure.structureMolecularWeight.max>70000.0</mvStructure.structureMolecularWeight.max>
  </orgPdbQuery>
 </queryRefinement>
 <queryRefinement>
  <queryRefinementLevel>5</queryRefinementLevel>
  <conjunctionType>and</conjunctionType>
  <orgPdbQuery>
    <queryType>org.pdb.query.simple.HomologueEntityReductionQuery</queryType>
    <description>Representative Structures at 90 Percent Sequence Identity</description>
    <queryId>null</queryId>
    <identityCutoff>90</identityCutoff>
  </orgPdbQuery>
 </queryRefinement>
</orgPdbCompositeQuery>

"""

### RUN QUERY AND RETURN RESULTS ##################################################################
req = urllib2.Request(url, data=queryText)

f = urllib2.urlopen(req)

for line in f:
    print line[0:4]
