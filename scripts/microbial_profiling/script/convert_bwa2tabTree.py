#!/usr/bin/env python3
import sys
from collections import defaultdict
from gi2lineage import *

# Load taxonomy
loadTaxonomy()

# Define major taxonomic levels
major_level = {
    'superkingdom': 10,
    'phylum': 20,
    'class': 30,
    'order': 40,
    'family': 50,
    'genus': 60,
    'species': 70
}

headers = []
taxa = defaultdict(int)

	#0 Accession
	#1 Length
	#2 GC%
	#3 Avg_fold
	#4 Fold_std
	#5 Base_Coverage%
	#6 Mapped_reads
	#7 Linear_length
	#8 ID


for line in sys.stdin:
    line = line.strip()
    if line.startswith("#") or not line:
        continue

    fields = line.split("\t")

    if fields[0] == "Accession":
        headers = fields[0:8]
        continue

    acc = fields[0]
    taxID = getTaxIDFromAcc(acc)

    if not taxID:
        print("[WARNING] Can't find Accession#{}: {}".format(acc, fields[8]), file=sys.stderr)
        continue

    ranks = []
    upper_rank_name = "root"
    for rank in sorted(major_level, key=major_level.get):
        name = acc2rank(acc, rank)
        name = name or "no rank - {}".format(upper_rank_name)
        name = name.replace("(", "").replace(")", "").replace(",", "").replace(":", "")
        ranks.append(name)
        upper_rank_name = name

    out = "\t".join(ranks)
    taxa[out] += int(fields[6])

for lineage in sorted(taxa, key=taxa.get, reverse=True):
    print("{}\t{}".format(taxa[lineage], lineage))
