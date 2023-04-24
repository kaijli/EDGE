#!/usr/bin/env python3

import sys
from gi2lineage import loadTaxonomy, getTaxRank, getTaxName, getTaxParent

min_count = int(sys.argv[1]) if len(sys.argv) > 1 else 0

# Load taxonomy
loadTaxonomy()

major_level = {
    'root': 0,
    'superkingdom': 10,
    'phylum': 20,
    'class': 30,
    'order': 40,
    'family': 50,
    'genus': 60,
    'species': 70,
    'strain': 80,
    'misc': 90
}

taxID = 0
taxa = {}

for line in sys.stdin:
    line = line.rstrip()
    if line.startswith('#') or line == '':
        continue

    if line.startswith('taxaID: '):
        taxID = int(line.split(': ')[1])
        continue

    if taxID > 1 and line.startswith('\tread count: '):
        read_count = int(line.split(': ')[1])

        if read_count < min_count:
            taxID = 0
            continue

        if taxID:
            rank = getTaxRank(taxID)

            # Deal with non-NCBI taxid
            if not rank:
                taxa['misc']['non-ncbi taxid']['READ_COUNT'] = taxa['misc']['non-ncbi taxid'].get('READ_COUNT', 0) + read_count
                taxa['misc']['non-ncbi taxid']['ROLLUP'] = taxa['misc']['non-ncbi taxid'].get('ROLLUP', 0) + read_count
                continue

            rank = 'strain' if rank == 'no rank'
            name = getTaxName(taxID)

            p_id = taxID
            while rank not in major_level:
                p_id = getTaxParent(p_id)
                rank = getTaxRank(p_id)
                name = getTaxName(p_id)

            taxa[rank][name]['READ_COUNT'] = taxa[rank][name].get('READ_COUNT', 0) + read_count

            while taxID:
                rank = getTaxRank(taxID)
                rank = 'strain' if rank == 'no rank'
                name = getTaxName(taxID)
                if name == 'root':
                    break

                if rank in major_level:
                    taxa[rank][name]['ROLLUP'] = taxa[rank][name].get('ROLLUP', 0) + read_count
                    taxa[rank][name]['TAXID'] = taxID

                taxID = getTaxParent(taxID)

        taxID = 0

print("LEVEL\tTAXA\tROLLUP\tASSIGNED\tTAXID")
for rank in sorted(major_level, key=major_level.get):
    for name in sorted(taxa[rank], key=lambda x: taxa[rank][x]['ROLLUP'], reverse=True):
        result = '\t'.join(taxa[rank][name]['RESULT']) if 'RESULT' in taxa[rank][name] else ''
        print(f"{rank}\t{name}\t{taxa[rank][name].get('ROLLUP', 0)}\t{taxa[rank][name].get('READ_COUNT', '')}\t{taxa[rank][name].get('TAXID', '')}")
