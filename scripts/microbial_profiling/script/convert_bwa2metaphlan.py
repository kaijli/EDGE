import sys
import Storable

sqd_in = sys.argv[1] if len(sys.argv) > 1 else None

dn_lvl = {
    'k': 'p',
    'p': 'c',
    'c': 'o',
    'o': 'f',
    'f': 'g',
    'g': 's'
}


def travel(ref, idx, string):
    for rank in ref[idx]:
        count = ref[idx][rank]['COUNT']
        if count > 0 and idx in dn_lvl:
            next_idx = dn_lvl[idx]
            next_ref = ref[idx][rank]
            travel(next_ref, next_idx, f"{string}|{rank}")
        out = f"{string}|{rank}"
        out = out.lstrip('|')
        #print(out, "\t", (count/tol)*100, "\n") if count > 0 and type == 'pcnt' else None
        print(out, "\t", count, "\n") if count > 0

def getLineage(gi):
    taxid = None
    strain = None
    lineage = []

    # Find lineage in species tree
    for tid in speciesTree:
        taxid = tid
        if 'GI' in speciesTree[taxid]:
            for rep in speciesTree[taxid]['GI']:
                for db in speciesTree[taxid]['GI'][rep]:
                    if gi in speciesTree[taxid]['GI'][rep][db]:
                        strain = speciesTree[taxid]['GI'][rep][db][gi]
                        break

    # If found, add lineage into an array
    if strain is not None:
        ranks = ['SK', 'P', 'C', 'O', 'F', 'G', 'S']
        for rank in ranks:
            name = speciesTree[taxid][rank]
            if rank == 'S':
                for species in speciesTree[taxid][rank]:
                    if speciesTree[taxid][rank][species] == 'scientific name':
                        name = species
            rank = rank.replace('SK', 'K')
            name = name.replace(' ', '_')
            lineage.append(f"{rank}__{name.lower()}")

    out = '|'.join(lineage)
    return out


print("Loading species tree...", file=sys.stderr)
speciesTree = Storable.retrieve("/lato/traceyf/db/custom/2013-03/speciesTreeGI.dmp")
# print("Loading genome tree...")
# genomeVitals = Storable.retrieve("/lato/traceyf/db/custom/2013-03/genomeVitals.dmp")

count = 1

print("Parsing node2phylo output file...", file=sys.stderr)
with open(sqd_in, 'r') as NODES:
    for line in NODES:
        line = line.strip()
        if line.startswith('GI'):
            continue
        temp = line.split("\t")
        name = temp[0]
        phylo = getLineage(temp[0])
        ranks = phylo.split('|')

        # initial ref
        ref = phylo_all
        for rank in ranks:
            idx = rank[0] if rank[0].isalpha() else None
            if not idx:
                raise Exception(f"\nERROR: Can't parse taxon - {name}\n")
            node_phylo[name][idx] = rank
            ref[idx][rank]['COUNT'] = 0

            if idx == 's':
                break

            # create an unclassified down-level
            next_idx = dn_lvl[idx]
            next_rank = rank.replace(f'{idx}__', f'{next_idx}__') + "_unclassified"
            ref[idx][rank][next_idx][next_rank]['COUNT'] = 0

            # reference
            ref = ref[idx][rank][next_idx]

        print(f"{count}\r", end='', file=sys.stderr)
        count += 1
print("Done", file=sys.stderr)

print("Parsing bwa output to tree...", file=sys.stderr)
with open(sqd_in, 'r') as MCV:
    for line in MCV:
        line = line.strip()
        if line.startswith('#') or 'ID' in line or not line:
            continue
        temp = line.split("\t")
        nid = temp[0]
        assigned_num = temp[6]

        ref = phylo_all
        for idx in ('k', 'p', 'c', 'o', 'f', 'g', 's'):
            if not node_phylo[nid][idx]:
                break

            rank = node_phylo[nid][idx]
            ref[idx][rank]['COUNT'] += assigned_num

            next_idx = dn_lvl[idx]
            if next_idx and not node_phylo[nid][next_idx]:
                next_rank = rank.replace(f'{idx}__', f'{next_idx}__') + "_unclassified"
                ref[idx][rank][next_idx][next_rank]['COUNT'] += assigned_num

            ref = ref[idx][rank]

print("Done", file=sys.stderr)
print("Output results...", file=sys.stderr)

tol = 0
for idx in phylo_all['k']:
    tol += phylo_all['k'][idx]['COUNT']



travel(phylo_all, 'k', "")
print("Done", file=sys.stderr)
print("Total reads assigned:", tol, file=sys.stderr)
