#!/usr/bin/env python3
import sys

proj = sys.argv[1] if len(sys.argv) > 1 else None
rank = sys.argv[2] if len(sys.argv) > 2 else None
contig = {}

for line in sys.stdin:
    temp = line.strip().split("\t")
    if rank and temp[1] != rank:
        continue
    if temp[0] in contig and temp[1] in contig[temp[0]]:
        continue
    contig[temp[0]] = {temp[1]: 1}
    out = line
    if proj:
        out = proj + "\t" + out.lstrip("#")[1:] if proj else out
    print(out, end='')
