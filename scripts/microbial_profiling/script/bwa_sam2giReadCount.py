#!/usr/bin/env python

# This script is used to output read mapping results in 2 columns: READ_NAME and TAXONOMY.
# 
# USAGE:
#   
#   bwa_sam2read_taxa.py ([RANK]) ("preload") < [sam file] > [output]
#  
# The default RANK is species. That can be changed to any rank. The "preload" option will load the
# entire GI-taxonomy mapping table into memory. That can increase the speed if your input is a
# very large SAM file (>10G).
#
# EXAMPLE:
#   
#   samtools view sample1.bam | bwa_sam2classification.py > output.classification
#   cat sample2.sam | bwa_sam2classification.py genus preload > output.genus_classification
# 
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.
# 2014/02/21

import sys
from gi2lineage import getAccFromSeqID

list = {}
score_cutoff = int(sys.argv[1]) if len(sys.argv) > 1 else 0

for line in sys.stdin:
    line = line.strip()
    if line.startswith('@'):
        continue
    fields = line.split('\t')
    flag = int(fields[1])
    seqid = fields[2]
    as_value = None
    xs_value = None
    for field in fields[11:]:
        if field.startswith('AS:i:'):
            as_value = int(field.split(':')[-1])
        elif field.startswith('XS:i:'):
            xs_value = int(field.split(':')[-1])
    if flag & 4 == 0:
        if as_value is not None and (score_cutoff == 0 or as_value >= score_cutoff):
            acc = getAccFromSeqID(seqid)
            gi = acc
            list.setdefault(gi, {}).setdefault('MAPPED', 0)
            list[gi]['MAPPED'] += 1
            if xs_value is None or as_value > xs_value:
                list.setdefault(gi, {}).setdefault('UNIQUE', 0)
                list[gi]['UNIQUE'] += 1

for gi in list:
    list[gi].setdefault('UNIQUE', 0)
    print("{}\t{}\t{}".format(gi, list[gi]['MAPPED'], list[gi]['UNIQUE']))
