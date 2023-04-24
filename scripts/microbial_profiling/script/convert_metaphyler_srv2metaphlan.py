#!/usr/bin/env python
import sys

def travel(ref, idx, string):
    if idx in ref:
        for rank in ref[idx]:
            count = ref[idx][rank]['COUNT']
            if count > 0 and idx in dn_lvl:
                next_idx = dn_lvl[idx]
                next_ref = ref[idx][rank]
                travel(next_ref, next_idx, f"{string}|{rank}")
            out = f"{string}|{rank}"
            out = out.lstrip('|')
            # print(out, (count / tol) * 100) if count > 0
            print(out, count) if count > 0

node_phylo = {}
phylo_all = {}
cutoff = 0.9

if len(sys.argv) < 2:
    print("USAGE: {} <metaphyler_output> [<cutoff>]".format(sys.argv[0]), file=sys.stderr)
    sys.exit(1)

infile = sys.argv[1]
if len(sys.argv) >= 3:
    cutoff = float(sys.argv[2])

dn_lvl = {
    'k': 'p',
    'p': 'c',
    'c': 'o',
    'o': 'f',
    'f': 'g',
    'g': 's'
}

lvl = ['g', 'f', 'o', 'c', 'p']

print("Parsing metaphylerSRV output file...", file=sys.stderr)
with open(infile) as nodes:
    for line in nodes:
        line = line.strip()
        if line.startswith('@'):
            #"prefix".classify.tab
            #column 0: query sequence id
            #column 2: phylogenetic marker gene name
            #    column 2: best reference gene hit
            #    column 3: % similarity with best hit
            #    column 4: classification rule
            #column 5-9: taxonomic label at genus,family,order,class,phylum level

            #"prefix".genus|family|order|class|phylum.tab
            #    column 1: taxonomic clade name
            #    column 2: % relative abundances
            #    column 3: depth of coverage of genomes
            #    column 4: number of sequences binned to this clade
            #    column 5: similarity with reference genes (only available at the genus level)	
            continue
        temp = line.split('\t')
        taxas = temp[5:10]
        name = taxas[0]
        ranks = []
        ranks.append("k__Bacteria")
        for i in range(4, -1, -1):
            if taxas[i].startswith('{'):
                break
            taxas[i] = taxas[i].replace(' ({})'.format(temp[2]), '')
            taxas[i] = taxas[i].replace(' ', '_')
            ranks.append("{}__{}".format(lvl[i], taxas[i]))

        # initial ref
        ref = phylo_all
        for rank in ranks:
            idx = rank.split('__')[0]
            if idx == 's':
                break

            if not idx:
                raise Exception("ERROR: Can't parse taxon - {}".format(name))

            if name not in node_phylo:
                node_phylo[name] = {}

            node_phylo[name][idx] = rank
            if idx not in ref:
                ref[idx] = {}
            ref[idx][rank] = {'COUNT': 0}

            # create an unclassified down-level
            next_idx = dn_lvl[idx]
            next_rank = rank.replace('{}__'.format(idx), '{}__'.format(next_idx))
            next_rank += "_unclassified"
            ref[idx][rank][next_idx] = {next_rank: {'COUNT': 0}}

            # reference
            ref = ref[idx][rank]

print("Done", file=sys.stderr)

print("Parsing metaphyler classification file...", file=sys.stderr)
with open(infile) as mpyl:
    for line in mpyl:
        line = line.strip()
        if line.startswith('#') or not line:
            continue

        temp = line.split('\t')
        tid, score = temp[5], 1
        for tid_scr in temp[5:10]:
            if tid_scr.startswith('{'):
                continue
            nid, assigned_num = tid, score

            ref = phylo_all
            for idx in 'kpcofg':
                if nid not in node_phylo or idx not in node_phylo[nid]:
                    continue

                rank = node_phylo[nid][idx]
                ref[idx][rank]['COUNT'] += assigned_num

                next_idx = dn_lvl[idx]
                if next_idx not in node_phylo[nid]:
                    next_rank = rank.replace('{}__'.format(idx), '{}__'.format(next_idx))
                    next_rank += "_unclassified"
                    ref[idx][rank][next_idx][next_rank]['COUNT'] += assigned_num

                ref = ref[idx][rank]

            break



print("Done", file=sys.stderr)

print("Parsing metaphyler classification file...", file=sys.stderr)
infile = "input_file.txt" # Update with your input file name
phylo_all = {}
node_phylo = {}
dn_lvl = {'k': 'p', 'p': 'c', 'c': 'o', 'o': 'f', 'f': 'g'}

with open(infile, 'r') as MPYL:
    for line in MPYL:
        line = line.strip()
        if line.startswith('#') or not line:
            continue
        if '@' in line:
            continue
        temp = line.split('\t')
        tid, score = temp[5], 1
        for tid_scr in temp[5:10]:
            if not re.search(r'\{\w+\}', tid_scr):
                nid, assigned_num = tid, score
                ref = phylo_all
                for idx in 'kpcofg':
                    if nid not in node_phylo or idx not in node_phylo[nid]:
                        break
                    rank = node_phylo[nid][idx]
                    ref.setdefault(idx, {}).setdefault(rank, {'COUNT': 0})
                    ref[idx][rank]['COUNT'] += assigned_num
                    next_idx = dn_lvl[idx]
                    if next_idx not in node_phylo.get(nid, {}):
                        next_rank = rank.replace(f'{idx}__', f'{next_idx}__') + '_unclassified'
                        ref.setdefault(idx, {}).setdefault(rank, {}).setdefault(next_idx, {}).setdefault(next_rank, {'COUNT': 0})
                        ref[idx][rank][next_idx][next_rank]['COUNT'] += assigned_num
                    ref = ref[idx][rank]
                break
print("Done", file=sys.stderr)

print("Output results...", file=sys.stderr)
tol = sum(phylo_all['k'][idx]['COUNT'] for idx in phylo_all['k'])
travel(phylo_all, 'k', '')
print("Done", file=sys.stderr)
print(f"Total read assigned: {tol}")

