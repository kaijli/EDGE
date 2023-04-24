#!/usr/bin/env python3
#
# This script is used to output read mapping results in 2 columns: READ_NAME and TAXONOMY.
# 
# USAGE:
#   
#   bwa_sam2read_taxa.py ([RANK]) ("preload") < [sam file] > [output]
#  
# The default RANK is species. That can be changed to any rank. The "preload" option will load the
# entire GI-taxonomy mapping table into memory. That can increase the speed if your input is a
# very large SAM file (>10G). Note that "preload" option can only work properly with Python 3.5 or 
# above.
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
#

import sys
import getopt
from gi2lineage import loadTaxonomy, getAccFromSeqID, acc2rank

rank = "species"
preload = None

opts, args = getopt.getopt(sys.argv[1:], "", ["rank=", "preload="])
for opt, arg in opts:
    if opt == "--rank":
        rank = arg
    elif opt == "--preload":
        preload = arg

loadTaxonomy(preload)

for line in sys.stdin:
    if line.startswith('@SQ'):
        continue
    fields = line.strip().split('\t')

    # unmapped
    if int(fields[1]) & 4:
        print(fields[0], '\t')
    else: # mapped
        # gi = fields[2].split('|')[1]
        acc = getAccFromSeqID(fields[2])
        name = acc2rank(acc, rank)
        if name:
            print(fields[0], '\t', name)
        else:
            print(fields[0], '\t', 'no', rank, '(', fields[2], ')')
