#!/usr/bin/env python3
# contig_classifier.py
# ver 0.1
# 2013/09/25
#
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.

import sys
import getopt
import os
import time

opt = {}
res, args = getopt.getopt(sys.argv[1:], 'i:l:', ['input=', 'line=', 'use_upper_name', 'help'])
for option, value in res:
    if option in ('-i', '--input'):
        opt['input'] = value
    elif option in ('-l', '--line'):
        opt['line'] = value
    elif option == '--use_upper_name':
        opt['use_upper_name'] = True
    elif option == '--help':
        opt['help'] = True

if 'help' in opt or 'input' not in opt or not os.path.exists(opt['input']):
    print("Usage: {} [OPTIONS] --input <FILE>".format(sys.argv[0]))
    print("Options:")
    print("  --input | -i <STRING>  the sequence of contigs in FASTA format")
    print("  --line | -l <NUM>      split over lines")
    print("  --use_upper_name       split over contig/assembly name")
    print("                          (E.g.: contig100_1 -> contig100)")
    print("  --help                 display this help")
    sys.exit()

INPUT = opt['input']
LINE = int(opt['line']) if 'line' in opt else 20000

time_now = time.time()
part = 0
flag = 0
cnt = 0
ctg = ""
fh = None

with open(INPUT, 'r') as input_file:
    for line in input_file:
        temp = line.strip().split("\t")
        name = temp[0]
        if 'use_upper_name' in opt:
            name = name.split("_")[0]
        if name != ctg and cnt >= LINE:
            cnt = 0
            if fh:
                fh.close()
            part += 1
            fh = open("{}.part-{}".format(INPUT, part), 'w')
        ctg = name
        if fh:
            fh.write(line)
        cnt += 1

if fh:
    fh.close()
