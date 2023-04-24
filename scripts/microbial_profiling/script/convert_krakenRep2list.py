#!/usr/bin/env python
import os
import sys
import time

# Load KronaTools library
os.system('ktGetLibPath')
import KronaTools

# Set output buffering
count = 0
taxa = None
period = None
time = time.time()
period = KronaTools.timeInterval(time)
KronaTools.loadTaxonomy()
period = KronaTools.timeInterval(time)

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
    'replicon': 95
}

level = {
    'root': 0,
    'superkingdom': 14,
    'kingdom': 15,
    'subkingdom': 16,
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


def timeInterval(now):
    now = int(time.time() - now)
    return '{:02d}:{:02d}:{:02d}'.format(int(now / 3600), int((now % 3600) / 60), int(now % 60))


for line in sys.stdin:
    count += 1
    if count % 5000 == 0:
        period = timeInterval(time)
        number = format(count, ',d')
        print("[{}] {} taxonomy processed.".format(period, number), end='\r', file=sys.stderr)

    fields = line.strip().split('\t')
    if fields[1] == '0':
        continue

    taxID = fields[4]
    name = fields[5].lstrip()

    rank = getTaxRank(taxID) or 'no rank'
    depth = getTaxDepth(taxID) or '10'

    rank = "unclassified" if taxID == '0' else rank
    rank = "root" if taxID == '1' else rank
    rank = "strain" if rank == "no rank" and int(depth) > 6 else rank

    # skip extra strain levels
    # no rank -- Escherichia coli O111:H-
    # no rank --  Escherichia coli O111:H- str. 11128
    if rank == "strain" and fields[2] == '0':
        continue

    taxa.setdefault(rank, {}).setdefault(name, {})['READ_COUNT'] = fields[2]
    taxa.setdefault(rank, {}).setdefault(name, {})['ROLLUP'] = fields[1]
    taxa.setdefault(rank, {}).setdefault(name, {})['TAXID'] = taxID

period = timeInterval(time)
number = format(count, ',d')
print("[{}] {} sequences processed.".format(period, number), file=sys.stderr)

print("LEVEL\tTAXA\tROLLUP\tASSIGNED\tTAXID")
for rank in sorted(taxa.keys(), key=lambda x: level[x]):
    for name in sorted(taxa[rank].keys(), key=lambda x: taxa[rank][x]['ROLLUP'], reverse=True):
        print("{}\t{}\t{}\t{}\t{}".format(
            rank,
            name,
            taxa[rank][name]['ROLLUP'],
            taxa[rank][name]['READ_COUNT'] if 'READ_COUNT' in taxa[rank][name] else "",
            taxa[rank][name]['TAXID']
        ))
