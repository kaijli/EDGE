# gi2lineage: This is a Python module to convert gi/NCBI taxonomy ID
# to taxonomy id/name/lineage. Some codes are taken from the module
# from Krona (Brian Ondov, Nicholas Bergman, and Adam Phillippy
# , Battelle National Biodefense Institute (BNBI))
#
# Po-E (Paul) Li
# B-11, LANL
# 04/2013 v0.1
#
# Changes log:
# 04/07/2017 - Add getAccFromSeqID, getTaxIDFromAcc acc2taxID acc2name acc2rank acc2lineage
# 04/04/2014 - Add taxid2lineage method
#
# 10/29/2013 - fixed bug retrieving GI>536870909 (2GB limitation)
#              NOTE: THIS PACKAGE REQUIRES PERL v5.16 or ABOVE.
#
# 10/28/2013 - able to map old GI to updated GI by providing a mapping
#              table using loadUpdatedGI().
#              (eg: /users/218817/scratch/opt/src/krona/taxonomy/gi_updated.tab).

import os
import csv

#############
# Exported #
############

# List of functions to be exported from the module
__all__ = [
    'acc2taxID',
    'acc2name',
    'acc2rank',
    'acc2lineage',
    'gi2taxID',
    'gi2name',
    'gi2rank',
    'gi2rank_taxid',
    'taxid2rank',
    'taxid2rank_taxid',
    'taxid2name',
    'taxid2lineage',
    'gi2lineage',
    'getTaxRank',
    'getTaxName',
    'getTaxDepth',
    'getTaxIDFromGI',
    'getTaxIDFromAcc',
    'getAccFromSeqID',
    'getTaxParent',
    'loadTaxonomy',
    'loadUpdatedGI'
]

####################
# Global constants #
####################

libPath = os.popen('ktGetLibPath').read().strip()
taxonomyDir = os.path.join(libPath, '..', 'taxonomy')
fileTaxByAcc = 'all.accession2taxid.sorted'
DEBUG = 0

#################
# Lookup tables #
#################
taxDepths = []
taxParents = []
taxRanks = []
taxNames = []
updatedGI = {}
taxIDByAcc = {}
taxIDByGI = {}
taxIDByGIStr = ''
invalidAccs = {}
missingAccs = {}

#################################################################
# Function: acc2taxID
# Description: Input an accession number and return its taxonomy ID.
# Arguments: acc - accession number
# Returns: taxID - taxonomy ID
#################################################################
def acc2taxID(acc):
    checkTaxonomy()
    return getTaxIDFromAcc(acc)

#################################################################
# Function: acc2name
# Description: Input an accession number and return its taxonomy name.
# Arguments: acc - accession number
# Returns: name - taxonomy name
#################################################################
def acc2name(acc):
    checkTaxonomy()
    taxID = getTaxIDFromAcc(acc)
    name = getTaxName(taxID)
    return name

#################################################################
# Function: acc2rank
# Description: Convert accession number to a specific taxonomy rank.
# Arguments: acc - accession number, r - taxonomy rank
# Returns: name - taxonomy name
# Example:
#   acc2rank("NC_000964","genus")
#   Returns: Bacillus
#################################################################
def acc2rank(acc, r):
    checkTaxonomy()
    rank = ""
    taxID = getTaxIDFromAcc(acc)
    rank = getTaxRank(taxID)
    name = getTaxName(taxID)

    if taxID == 1:
        return "root"
    if r == "root":
        return "root"
    if r.lower() == "strain":
        return name

    while taxID:
        if rank.lower() == r.lower():
            return name

        if name == 'root':
            break

        taxID = getTaxParent(taxID)
        rank = getTaxRank(taxID)
        name = getTaxName(taxID)

    return ""

#################################################################
# Function: gi2taxID
# Description: Input a GI number and return its taxonomy ID.
# Arguments: gi (GI number)
# Return: Taxonomy ID
#################################################################
def gi2taxID(gi):
    # Call checkTaxonomy() function
    checkTaxonomy()
    # Call getTaxIDFromGI() function with gi as argument and return the result
    return getTaxIDFromGI(gi)

#################################################################
# Function: gi2name
# Description: Input a GI number and return its taxonomy name.
# Arguments: gi (GI number)
# Return: Taxonomy name
#################################################################
def gi2name(gi):
    # Call checkTaxonomy() function
    checkTaxonomy()
    # Call getTaxIDFromGI() function with gi as argument to get taxID
    taxID = getTaxIDFromGI(gi)
    # Call getTaxName() function with taxID as argument and return the result
    name = getTaxName(taxID)
    return name

#################################################################
# Function: taxid2name
# Description: Input a taxonomy ID and return its taxonomy name.
# Arguments: taxID (taxonomy ID)
# Return: Taxonomy name
#################################################################
def taxid2name(taxID):
    # Call checkTaxonomy() function
    checkTaxonomy()
    # Call getTaxName() function with taxID as argument and return the result
    name = getTaxName(taxID)
    return name

#################################################################
# Function: gi2rank
# Description: Convert GI to a specific rank.
# Arguments: gi (GI number), r (desired rank)
# Return: Taxonomy name
# Example:
#      gi2rank("47552137","genus")
#      Returns: Bacillus
#################################################################
def gi2rank(gi, r):
    # Call checkTaxonomy() function
    checkTaxonomy()
    # Initialize variables
    rank = ""
    name = ""
    taxID = getTaxIDFromGI(gi)
    rank = getTaxRank(taxID)
    name = getTaxName(taxID)
    
    # Return "root" if taxID is 1 or if r is "root"
    if taxID == 1 or r.lower() == "root":
        return "root"
    # Return name if r is "strain"
    elif r.lower() == "strain":
        return name
    
    # Loop until taxID becomes 0 or name becomes "root"
    while taxID:
        # Return name if rank matches r (case insensitive)
        if rank.lower() == r.lower():
            return name
        # Break if name is "root"
        if name == "root":
            break
        # Update taxID, rank, and name to parent values
        taxID = getTaxParent(taxID)
        rank = getTaxRank(taxID)
        name = getTaxName(taxID)
    
    return ""

#################################################################
# Function: taxid2rank
# Description: Convert taxonomy ID to a specific rank.
# Arguments: id (taxonomy ID), r (desired rank)
# Return: Taxonomy name
#################################################################
def taxid2rank(id, r):
    # Call checkTaxonomy() function
    checkTaxonomy()
    # Return "root" if id is 1 or if r is "root"
    if id == 1 or r.lower() == "root":
        return "root"
    # Initialize variables
    rank = ""
    name = ""
    taxID = id
    rank = getTaxRank(taxID) or 'no rank'
    name = getTaxName(taxID) or ''
    
    # Return name if r is "strain"
    if r.lower() == "strain":
        return name
    
    # Loop until taxID becomes 0 or name becomes "root"
    while taxID:
        # Return name if rank matches r (case insensitive)
        if rank.lower() == r.lower():
            return name
        # Break if name is "root"
        if name == "root":
            break
        # Update taxID, rank, and name to parent values
        taxID = getTaxParent(taxID)
        rank = getTaxRank(taxID)
        name = getTaxName(taxID)
    
    return ""

#################################################################
# Function: gi2rank_taxid
# Description: Convert GI to a specific taxonomic rank
# Arguments: gi - GI number
#            r - Target taxonomic rank
# Returns: Taxonomy ID
# Example: gi2rank_taxid("47552137", "genus") returns "Bacillus"
#################################################################
def gi2rank_taxid(gi, r):
    checkTaxonomy() # Call to checkTaxonomy() function (not defined in the code)

    # Initialize variables
    rank = ""
    taxID = getTaxIDFromGI(gi) # Call to getTaxIDFromGI() function (not defined in the code)
    rank = getTaxRank(taxID) # Call to getTaxRank() function (not defined in the code)
    name = getTaxName(taxID) # Call to getTaxName() function (not defined in the code)

    # Return taxID if rank is 'strain'
    if rank == 'strain':
        return taxID

    # Return 1 if rank is 'root'
    if rank == 'root':
        return 1
    
    # Loop until taxID becomes empty
    while taxID:
        # Return taxID if rank matches the target rank (case insensitive)
        if rank.lower() == r.lower():
            return taxID
        
        # Exit loop if name is 'root'
        if name == 'root':
            break
        
        # Update taxID, rank, and name using getTaxParent(), getTaxRank(), and getTaxName() functions respectively
        taxID = getTaxParent(taxID) # Call to getTaxParent() function (not defined in the code)
        rank = getTaxRank(taxID)
        name = getTaxName(taxID)
    
    # Return empty string if target rank not found
    return ""

#################################################################
# Function: taxid2rank_taxid
# Description: Convert Taxonomy ID to a specific taxonomic rank
# Arguments: id - Taxonomy ID
#            r - Target taxonomic rank
# Returns: Taxonomy ID
#################################################################
def taxid2rank_taxid(id, r):
    checkTaxonomy()  # Call the checkTaxonomy() function

    taxID = id  # Set the initial value of taxID to the input id
    rank = getTaxRank(taxID)  # Get the taxonomic rank of the given taxID
    name = getTaxName(taxID)  # Get the taxonomic name of the given taxID

    # Check if the given rank matches the rank of the taxID, or if the given rank is 'strain' and the rank of the taxID is 'no rank'
    if r == rank or (r == 'strain' and rank == 'no rank'):
        return taxID  # Return the taxID

    if r == "root":  # Check if the given rank is 'root'
        return 1  # Return 1

    while taxID:
        if rank.upper() == r.upper():  # Check if the rank of the taxID matches the given rank (case-insensitive)
            return taxID  # Return the taxID
        if name == 'root':  # Check if the taxonomic name is 'root'
            break  # Exit the loop

        taxID = getTaxParent(taxID)  # Get the parent taxID of the current taxID
        rank = getTaxRank(taxID)  # Get the taxonomic rank of the current taxID
        name = getTaxName(taxID)  # Get the taxonomic name of the current taxID

    return ""  # Return an empty string if no matching taxID is found

#################################################################
# function: gi2lineage(<GI>[,<PRINT_ALL_RANK>,<PRINT_STRAIN>,<REMOVE_SPACE>])
# input:    <GI> - GI number
#           <PRINT_ALL_RANK> (1|0) - (optional, default is 1)
#               1: display 7 major ranks even they are not available (display as no_rank).
#               0: display available ranks only.
#           <PRINT_STRAIN> (1|0) - (optional, default is 0)
#               1: display guessed strain name (n).
#               0: do not display strain name.
#           <REMOVE_SPACE> (1|0) - (optional, default is 1)
#               1: remove space from organism name to "_"
#               0: do not remove space from organism name
# output:   The full lineage string is joined of the ranks with "|" separator. Each rank represents by its
#           initial as an abbreviation following by two underscores ("__") and the taxonomy.
#
# example:  gi2lineage("47552137")
#           returns: k__Bacteria|p__Firmicutes|c__Bacilli|o__Bacillales|f__Bacillaceae|g__Bacillus|s__Bacillus_anthracis
#################################################################

def taxid2lineage(id, print_all_rank=1, print_strain=0, replace_space2underscore=1):
    id2lineage(id, print_all_rank, print_strain, replace_space2underscore, "taxid")

def gi2lineage(id, print_all_rank=1, print_strain=0, replace_space2underscore=1):
    id2lineage(id, print_all_rank, print_strain, replace_space2underscore, "gi")

def acc2lineage(id, print_all_rank=1, print_strain=0, replace_space2underscore=1):
    id2lineage(id, print_all_rank, print_strain, replace_space2underscore, "acc")

def id2lineage(id, print_all_rank=None, print_strain=None, replace_space2underscore=None, input_type=None):
    # Check if taxonomy data is available
    checkTaxonomy()

    # Set default values if not defined
    print_all_rank = 1 if print_all_rank is None else print_all_rank
    print_strain = 0 if print_strain is None else print_strain
    replace_space2underscore = 1 if replace_space2underscore is None else replace_space2underscore

    # Define major taxonomic levels and their corresponding codes
    major_level = {
        'superkingdom': 'k',
        'phylum': 'p',
        'class': 'c',
        'order': 'o',
        'family': 'f',
        'genus': 'g',
        'species': 's'
    }

    # Define level codes and initialize with empty strings
    level = {
        'k': '',
        'p': '',
        'c': '',
        'o': '',
        'f': '',
        'g': '',
        's': ''
    }

    # Get taxonomic ID based on input type
    if input_type == 'taxid':
        taxID = id
    elif input_type == 'gi':
        taxID = getTaxIDFromGI(id)
    elif input_type == 'acc':
        taxID = getTaxIDFromAcc(id)

    # Get initial taxonomic rank and name
    rank = getTaxRank(taxID)
    name = getTaxName(taxID)
    str_name = name
    if replace_space2underscore:
        str_name = str_name.replace(' ', '_')

    # Traverse up the taxonomic tree
    while taxID:
        if rank in major_level:
            if replace_space2underscore:
                name = name.replace(' ', '_')
            level[major_level[rank]] = name
        if name == 'root':
            break

        taxID = getTaxParent(taxID)
        rank = getTaxRank(taxID)
        name = getTaxName(taxID)

    last = str_name
    rank_list = []
    # Generate the lineage string for each taxonomic level
    for lvl in ['s', 'g', 'f', 'o', 'c', 'p', 'k']:
        if not print_all_rank and not level[lvl]:
            continue
        level[lvl] = f"{last} - no_{lvl}_rank" if not level[lvl] else level[lvl]
        last = level[lvl]
        rank_list.append(f"{lvl}__{level[lvl]}")

    # Reverse the rank list to get the correct order
    rank_list.reverse()

    # Add strain information if specified
    if print_strain:
        rank_list.append(f"n__{str_name}")

    # Join the rank list with '|' separator
    return '|'.join(rank_list)

def getTaxDepth(taxID):
    return taxDepths[taxID]

def getTaxName(taxID):
    return taxNames[taxID]

def getTaxParent(taxID):
    while taxID > 1 and taxRanks[taxID] == 'no rank':
        taxID = taxParents[taxID] or 1
    return taxID

def getTaxRank(taxID):
    return taxRanks[taxID]

def getTaxIDFromGI(gi):
    gi = updatedGI[gi] if gi in updatedGI else gi
    taxID = None

    if gi in taxIDByGI:
        return taxIDByGI[gi]
    elif taxIDByGIStr is not None:
        pos = gi * 4
        try:
            # the last one x2147483636L (2GB limitation)
            # taxID = unpack("x${pos}L", taxIDByGIStr)
            data = taxIDByGIStr[pos:pos+4]
            taxID = struct.unpack("L", data)[0]
        except:
            raise Exception("ERROR: GI to TaxID table is out of date.\n")
    elif gi not in taxIDByGI:
        with open(f"{taxonomyDir}/gi_taxid.dat", "rb") as gi_file:
            gi_file.seek(gi * 4)
            data = gi_file.read(4)
            taxID = struct.unpack("L", data)[0]
    # make sure taxID is in database
    if taxID is not None and taxRanks[taxID] is not None:
        taxIDByGI[gi] = taxID

    return taxIDByGI[gi]

def getTaxIDFromAcc(acc):
    if acc.isdigit():
        return acc

    acc = acc.split('.')[0]
    
    if acc in taxIDByAcc:
        return taxIDByAcc[acc]
    
    size = os.path.getsize(f"{taxonomyDir}/{fileTaxByAcc}")
    accCur = ""
    taxID = None
    
    with open(f"{taxonomyDir}/{fileTaxByAcc}", "r") as ACC:
        if not ACC:
            print("ERROR: Sorted accession to taxID list not found. Was updateAccessions.sh run?")
            exit(1)
        
        min_pos = 0
        max_pos = size
        while acc != accCur and min_pos < max_pos:
            pos_new = int((min_pos + max_pos) / 2)
            
            ACC.seek(pos_new, 0)
            
            if pos_new != min_pos:
                ACC.readline()  # eat up to newline
            
            line = ACC.readline()
            
            acc_new, taxID = line.split('\t')
            
            if acc > acc_new and accCur != acc_new and acc_new:
                if acc_new:
                    pos_new = ACC.tell()
                
                min_pos = pos_new
                
                if min_pos >= max_pos:
                    max_pos = min_pos + 1
            else:
                max_pos = pos_new
            
            accCur = acc_new
        
        ACC.close()
    
    if accCur != acc:
        missingAccs[acc] = 1
        taxID = 0
    
    taxIDByAcc[acc] = taxID
    
    return taxIDByAcc[acc]

def getAccFromSeqID(seqID):
    acc = seqID.split()[0]
    
    if '|' in acc:
        acc = acc.split('|')[3]
    
    if not acc.isdigit() and not re.match(r"^[A-Z]+_?[A-Z]*\d+(\.\d+)?$", acc):
        invalidAccs[acc] = 1
        return None
    
    return acc

def load_updated_gi(gi_updated_file):
    updated_gi = {}
    with open(gi_updated_file, 'r') as info:
        for line in info:
            line = line.strip()
            old_gi, new_gi = line.split('\t')
            updated_gi[old_gi] = new_gi
    return updated_gi

def load_taxonomy(gi_tax_file=None):
    taxParents = []
    taxDepths = []
    taxRanks = []
    taxNames = []
    taxonomy_file = f"{taxonomyDir}/taxonomy.tab" 
    with open(taxonomy_file, 'r') as info:
        for line in info:
            line = line.strip()
            id, depth, parent, rank, name = line.split('\t')
            taxParents.append(parent)
            taxDepths.append(depth)
            taxRanks.append(rank)
            taxNames.append(name)
    print("Done parsing taxonomy.tab")  # Replace with appropriate debug statement

    if taxParents[2] == 1:
        print("Local taxonomy database is out of date. Update using updateTaxonomy.sh.")  # Replace with appropriate error message

    if gi_tax_file and os.path.exists(gi_tax_file):
        with open(gi_tax_file, 'rb') as f:
            tax_id_by_gi = pickle.load(f)
    elif gi_tax_file and not os.path.exists(gi_tax_file) and 'load' in gi_tax_file.lower():
        gi_tax_data_file = f"{taxonomyDir}/gi_taxid.dat" 
        with open(gi_tax_data_file, 'rb') as f:
            tax_id_by_gi = pickle.load(f)
    return taxParents, taxDepths, taxRanks, taxNames, tax_id_by_gi

def tax_contains(parent, child):
    # Determines if parent is an ancestor of (or equal to) child

    depth_parent = taxDepths[parent]

    while taxDepths[child] > taxDepths[parent]:
        child = taxParents[child]

    return parent == child


################
# Not exported #
################
def kt_die(error):
    print('[ ERROR ]', error)
    exit(1)


def check_taxonomy():
    if not taxParents:
        raise Exception('Taxonomy not loaded. "load_taxonomy()" must be called first.')


def taxon_link(tax_id):
    return tax_id


# Usage of the functions
# load_taxonomy() # Call this function to load the taxonomy data into the global variables
# tax_contains(parent, child) # Call this function to check if parent is an ancestor of child
# kt_die(error) # Call this function to print an error message and exit with error code 1
# check_taxonomy() # Call this function to check if the taxonomy data has been loaded
# taxon_link(tax_id) # Call this function to get a link for a given taxon ID
