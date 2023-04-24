#!/usr/bin/env python3

import sys
import re
from gi2lineage import loadTaxonomy, getTaxRank, getTaxDepth, taxid2rank

col = None

major_level = {
    'root': 0,
    'superkingdom': 15,
    'phylum': 25,
    'class': 35,
    'order': 45,
    'family': 55,
    'genus': 65,
    'species': 75,
    'strain': 85,
    # 'replicon': 95
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

for line in sys.stdin:
    line = line.rstrip()

    if re.search(r'^LEVEL', line) and not col:
        col = 0
        header_titles = line.split('\t')
        for header_line in header_titles:
            if re.search(r'^(TAXID|TAXAID|TAXOID|TAXA_ID)$', header_line):
                break
            col += 1
        continue

    fields = line.split('\t')

    taxID = fields[col]
    name = fields[1].lstrip()

    rank = getTaxRank(taxID)
    depth = getTaxDepth(taxID)

    if taxID == 0:
        rank = 'unclassified'
    if taxID == 1:
        rank = 'root'
    if rank == 'no rank' and depth > 6:
        rank = 'strain'

    # transfer minor rank to major rank
    close_major_rank = rank
    if rank not in major_level:
        for temp in sorted(major_level, key=major_level.get):
            if level[rank] > major_level[temp]:
                close_major_rank = temp
        pass

    lineage = []
    for rank in sorted(major_level, key=major_level.get):
        rname = taxid2rank(taxID, rank)
        if rank == 'strain':
            rname = name
        if not rname:
            rname = 'no rank'
        lineage.append(rname)
        if rank == close_major_rank:
            break

    lineage_idx = '\t'.join(lineage)
    lineage_idx = re.sub(r'([\(\),:])', '', lineage_idx)
    taxa[lineage_idx] = taxa.get(lineage_idx, 0) + int(fields[2])

for lineage_idx in sorted(taxa.keys()):
    score = taxa[lineage_idx]
    print(f"{score}\t{lineage_idx}")
