#!/usr/bin/env python3

import sys
from gi2lineage import *

loadTaxonomy()

cutoff = int(sys.argv[1]) if len(sys.argv) > 1 else 0
major_level = {
    'superkingdom': 10,
    'phylum': 20,
    'class': 30,
    'order': 40,
    'family': 50,
    'genus': 60,
    'species': 70
    # 'strain': 80,
    # 'replicon': 90
}

headers = []
taxa = {}

for line in sys.stdin:
    line = line.rstrip()
    if line.startswith("#") or not line:
        continue

    fields = line.split("\t")

    acc = fields[0]
    taxID = getTaxIDFromAcc(acc)

    if not taxID:
        print("[WARNING] Can't find Accession#{}".format(acc), file=sys.stderr)
        continue

    ranks = []
    upper_rank_name = "root"
    for rank in sorted(major_level, key=major_level.get):
        name = acc2rank(acc, rank)
        name = "no rank - {}".format(upper_rank_name) if not name else name
        name = name.replace("(", "").replace(")", "").replace(",", "").replace(":", "")
        ranks.append(name)
        upper_rank_name = name

    out = "\t".join(ranks)

    if out not in taxa:
        taxa[out] = 0

    taxa[out] += int(fields[1])

for lineage in sorted(taxa):
    if taxa[lineage] >= cutoff:
        print("{}\t{}".format(taxa[lineage], lineage))
