#!/usr/bin/env python
import sys

taxa = {}

for line in sys.stdin:
    line = line.strip()
    if "@" in line:
        continue
    temp = line.split("\t")
    taxas = list(reversed(temp[5:10]))
    lineage = "\t".join(taxas)
    lineage = lineage.replace("{\w+}", "").replace("([\(\),:])", "")
    taxa[lineage] = taxa.get(lineage, 0) + 1

for lineage in sorted(taxa.keys(), key=lambda x: taxa[x], reverse=True):
    print(f"{taxa[lineage]}\t{lineage}")
