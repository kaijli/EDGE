#!/usr/bin/env python
import sys

tax = {}

for line in sys.stdin:
    line = line.strip()
    if "LEVEL" in line:
        continue
    if "|s__" not in line:
        continue
    phylo, count = line.split("\t")
    phylo = phylo.replace(r"^\w__", "\t").replace(r"\|\w__", "\t").replace("_", " ")
    tax[phylo] = count

for phylo in sorted(tax.keys()):
    print(f"{tax[phylo]}\t{phylo}")
