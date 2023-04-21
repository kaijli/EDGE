#!/usr/bin/env python
import sys
import getopt

fasta = ""
taxa = ""
csv = ""
rank = "strain"

# Get command line options
opts, args = getopt.getopt(sys.argv[1:], "fasta:taxa:csv:rank:help?")
for opt, arg in opts:
    if opt == "-fasta":
        fasta = arg
    elif opt == "-taxa":
        taxa = arg
    elif opt == "-csv":
        csv = arg
    elif opt == "-rank":
        rank = arg
    elif opt in ("-help", "-?"):
        print("Usage: python {} -fasta <contig> -csv <ctg_class.top.csv> -taxa NameOftaxa -rank TaxRank".format(sys.argv[0]))
        print("e.g: python {} -fasta contigs.fasta -csv output.ctg_class.top.csv -taxa 'Zaire ebolavirus' -rank strain".format(sys.argv[0]))
        print("-rank    phylum, class, order, family, genus, species, or strain (default: strain)")
        sys.exit()

if not (fasta and taxa and csv):
    print("Error: -fasta, -taxa, and -csv options are required.")
    print("Usage: python {} -fasta <contig> -csv <ctg_class.top.csv> -taxa NameOftaxa -rank TaxRank".format(sys.argv[0]))
    sys.exit()

hit_id = {}
with open(csv, 'r') as fh:
    for line in fh:
        line = line.rstrip()
        tmp = line.split('\t')
        if tmp[1] == rank and taxa in tmp[2]:
            hit_id[tmp[0]] = tmp[2]

with open(fasta, 'r') as cfh:
    for line in cfh:
        line = line.strip()
        if not line.startswith('>'):
            continue
        head, *seq = line.split('\n')
        id, desc = head[1:].split(None, 1)
        if id in hit_id:
            id = "{} {}{}".format(id, desc, hit_id[id])
            seq = ''.join(seq)
            seq = '\n'.join([seq[i:i+70] for i in range(0, len(seq), 70)])
            print(">{}\n{}".format(id, seq))
