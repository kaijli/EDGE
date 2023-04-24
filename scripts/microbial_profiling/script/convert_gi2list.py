#!/usr/bin/env python3

from gi2lineage import *
from collections import defaultdict

# load taxonomy
print ("Loading taxonomy...\n", file=sys.stderr)
loadTaxonomy()

major_level = {
    'superkingdom': 10,
    'phylum': 20,
    'class': 30,
    'order': 40,
    'family': 50,
    'genus': 60,
    'species': 70,
    'strain': 80,
    'replicon': 90
}

headers = []
taxa = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

while True:
    line = input().rstrip()
    if not line:
        break

    fields = line.split("\t")

    acc = fields[0]
    taxID = getTaxIDFromAcc(acc)

    print(f"[WARNING] Can't find Accession#{acc}") if not taxID

    rank = "replicon"
    name = acc

    taxa[rank][name]['READ_COUNT'] += int(fields[1])
    taxa[rank][name]['ROLLUP'] += int(fields[1])

    if taxID:
        #print("Processing GI:$gi TAXID:$taxID NAME:'$fields[8]'...\n", file=sys.stderr);
        rank = getTaxRank(taxID) or "no rank"
        rank = "strain" if rank == "no rank"
        name = getTaxName(taxID)

        p_id = taxID
        while rank not in major_level:
            p_id = getTaxParent(p_id)
            rank = getTaxRank(p_id)
            name = getTaxName(p_id)

        taxa[rank][name]['READ_COUNT'] += int(fields[1])

        while taxID:
            rank = getTaxRank(taxID)
            rank = "strain" if rank == "no rank"
            name = getTaxName(taxID)

            if name == 'root':
                break

            if rank in major_level:
                taxa[rank][name]['ROLLUP'] += int(fields[1])

            taxID = getTaxParent(taxID)

print("\n", file=sys.stderr)
print("LEVEL\tTAXA\tROLLUP\tASSIGNED")
for rank in sorted(major_level, key=major_level.get):
    for name in sorted(taxa[rank], key=lambda x: taxa[rank][x]['ROLLUP'], reverse=True):
        print(f"{rank}\t{name}\t{taxa[rank][name]['ROLLUP']}\t{taxa[rank][name]['READ_COUNT'] if 'READ_COUNT' in taxa[rank][name] else ''}")
