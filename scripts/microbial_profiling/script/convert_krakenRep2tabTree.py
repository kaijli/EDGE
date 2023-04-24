#!/usr/bin/env python3
from gi2lineage import *

# load taxonomy
loadTaxonomy()

major_level = {
    'root': 0,
    'superkingdom': 15,
    'phylum': 25,
    'class': 35,
    'order': 45,
    'family': 55,
    'genus': 65,
    'species': 75,
    'strain': 85
    #'replicon'     => 95
}

level = {
    'root': 0,
    'superkingdom': 15,
    'kingdom': 16,
    'subkingdom': 17,
    'superphylum': 24,
    'phylum': 25,
    'subphylum': 26,
    'infraclass': 33,
    'superclass': 34,
    'class': 35,
    'subclass': 36,
    'infraorder': 43,
    'superorder': 44,
    'order': 45,
    'suborder': 46,
    'parvorder': 47,
    'superfamily': 54,
    'family': 55,
    'subfamily': 56,
    'tribe': 58,
    'genus': 65,
    'subgenus': 66,
    'species group': 73,
    'species subgroup': 74,
    'species': 75,
    'subspecies': 76,
    'strain': 85,
    'no rank': 90,
    'unclassified': 100
}

taxa = {}
count = 0

#  23.85  1564844 1564844 U   0   unclassified
#  76.15  4997221 7839    -   1   root
#  76.03  4989140 105 -   131567    cellular organisms
#  75.47  4952327 9220    -   2       Bacteria
#  28.79  1889195 0   P   1297          Deinococcus-Thermus
#  28.79  1889195 661 C   188787          Deinococci
#  28.78  1888527 106 O   118964            Deinococcales
#  28.78  1888419 0   F   183710              Deinococcaceae
#  28.78  1888419 14154   G   1298                  Deinococcus
#  28.56  1873957 0   S   1299                    Deinococcus radiodurans
#  28.56  1873957 1873957 -   243230                    Deinococcus radiodurans R1
#   0.00  0   0   -   1408434                   Deinococcus radiodurans ATCC 13939
#   0.00  145 0   S   502394                  Deinococcus gobiensis
#   0.00  145 145 -   745776                    Deinococcus gobiensis I-0
#   0.00  64  0   S   55148                   Deinococcus proteolyticus
#   0.00  64  64  -   693977                    Deinococcus proteolyticus MRP
#   0.00  38  0   S   310783                  Deinococcus deserti


while True:
    try:
        line = input().rstrip('\n')
    except EOFError:
        break

    fields = line.split('\t')
    if fields[1] == '0' or fields[2] == '0':
        continue

    taxID = fields[4]
    name = fields[5].lstrip()

    rank = getTaxRank(taxID) or 'no rank'
    depth = getTaxDepth(taxID) or '10'

    if taxID == '0':
        rank = 'unclassified'
    if taxID == '1':
        rank = 'root'
    if rank == 'no rank' and int(depth) > 6:
        rank = 'strain'

    # transfer minor rank to major rank
    close_major_rank = rank
    if rank not in major_level:
        for temp in sorted(major_level.keys(), key=lambda x: major_level[x]):
            if level[rank] > major_level[temp]:
                break
            close_major_rank = temp

    lineage = []
    for rank in sorted(major_level.keys(), key=lambda x: major_level[x]):
        rname = taxa[taxID][rank] if taxID in taxa and rank in taxa[taxID] else ''
        lineage.append(rname)
    print("\t".join(lineage))
