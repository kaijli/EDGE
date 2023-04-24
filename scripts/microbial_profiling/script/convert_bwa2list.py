#!/usr/bin/python
import sys
import os
import re
import KronaTools

# load taxonomy
print("Loading taxonomy...")
loadTaxonomy()

# parse BLAST results
fileName = sys.argv[1]

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
taxa = {}

with open(fileName, 'r') as TAXA:
	for line in TAXA:
		line = line.rstrip()
		if re.search('^#', line) or not line:
			continue

		# 0 Accession
		# 1 Length
		# 2 GC%
		# 3 Avg_fold
		# 4 Fold_std
		# 5 Base_Coverage%
		# 6 Mapped_reads
		# 7 Linear_length
		# 8 ID
		fields = line.split('\t')

		if fields[0] == "GI":
			headers = fields[0:7]
			continue

		acc = fields[0]
		taxID = getTaxIDFromGI(acc)

		if not taxID:
			print("[WARNING] Can't find Accession#{}: {}".format(acc, fields[8]), file=sys.stderr)

		rank = "replicon"
		name = fields[8]

		if rank not in taxa:
			taxa[rank] = {}
		if name not in taxa[rank]:
			taxa[rank][name] = {}
		if 'READ_COUNT' not in taxa[rank][name]:
			taxa[rank][name]['READ_COUNT'] = 0
		taxa[rank][name]['READ_COUNT'] += int(fields[6])
		if 'ROLLUP' not in taxa[rank][name]:
			taxa[rank][name]['ROLLUP'] = 0
		taxa[rank][name]['ROLLUP'] += int(fields[6])
		if 'RESULT' not in taxa[rank][name]:
			taxa[rank][name]['RESULT'] = fields[0:7]

		if taxID:
			rank = getTaxRank(taxID)
			rank = "strain" if rank == "no rank"
			name = getTaxName(taxID)

			p_id = taxID
			while rank not in major_level:
				p_id = getTaxParent(p_id)
				rank = getTaxRank(p_id)
				name = getTaxName(p_id)
			if rank not in taxa:
				taxa[rank] = {}
			if name not in taxa[rank]:
				taxa[rank][name] = {}
			if 'READ_COUNT' not in taxa[rank][name]:
				taxa[rank][name]['READ_COUNT'] = 0
			taxa[rank][name]['READ_COUNT'] += int(fields[6])

			while taxID:
                rank = getTaxRank(taxID)
                rank = "strain" if rank == "no rank" else rank

                name = getTaxName(taxID)
                if name == 'root':
                    break

                if rank in major_level:
                    taxa[rank][name]['ROLLUP'] = 0 if 'ROLLUP' not in taxa[rank][name] else taxa[rank][name]['ROLLUP']
                    taxa[rank][name]['ROLLUP'] += int(fields[6])

                taxID = getTaxParent(taxID)

# Close TAXA file
TAXA.close()

# Print a new line to STDERR
print("\n", file=sys.stderr)

# Prepare the header for output
header = '\t'.join(headers)
print(f"LEVEL\tTAXA\tROLLUP\tASSIGNED\t{header}")

# Loop through taxa dictionary and print results
for rank in sorted(major_level, key=major_level.get):
    for name in sorted(taxa[rank], key=lambda x: taxa[rank][x]['ROLLUP'], reverse=True):
        result = '\t'.join(taxa[rank][name]['RESULT']) if 'RESULT' in taxa[rank][name] else ''
        print(f"{rank}\t{name}\t{taxa[rank][name]['ROLLUP']}\t{taxa[rank][name]['READ_COUNT'] if 'READ_COUNT' in taxa[rank][name] else ''}\t{result}")

