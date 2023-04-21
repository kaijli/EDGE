#!/usr/bin/env python3

import argparse
from os.path import exists
from gi2lineage import loadTaxonomy
import sys
import re

def usage():
    print("""Usage: python gi2lineage.py -i input -t type -p orig_seq

Options:
  -i, --input=input   input file (required)
  -t, --type=type     input format: last, blast, sam (default: guess)
  -p, --orig_seq=file original query sequence file
  -h, --help          show this help message and exit""").format(sys.argv[0])
  sys.exit()

parser = argparse.ArgumentParser(description="Process GI to lineage information")
parser.add_argument("-i", "--input", dest="input", required=True,
                    help="input file")
parser.add_argument("-t", "--type", dest="type", default="guess",
                    help="input format: last, blast, sam (default: guess)")
parser.add_argument("-p", "--orig_seq", dest="orig_seq",
                    help="original query sequence file")
parser.add_argument("-?", "--help", dest="help", action='store_true',
                    help="show this help message and exit")

args = parser.parse_args()

if args.help or not exists(args.input):
    usage()
    exit()

print("Input file: " + args.input)

# Preload taxonomy data
# loadTaxonomy("preload")
loadTaxonomy()
print("Done loading taxonomy data.")

input_file = args.input
input_format = args.type
file = args.orig_seq

input_format = "guess" if not input_format else input_format

#############################################################################################################
## Support 3 alignment format:
#############################################################################################################
#
#  # LAST tab output    # BLAST M8       SAM Field  Description
#  0  score             0  qseqid        0   QNAME  Query (pair) NAME
#  1  name1 (ref)       1  sseqid        1   FLAG   bitwise FLAG
#  2  start1            2  pident        2   RNAME  Reference sequence NAME
#  3  alnSize1          3  length        3   POS    1-based leftmost POSition/coordinate of clipped sequence
#  4  strand1           4  mismatch      4   MAPQ   MAPping Quality (Phred-scaled)
#  5  seqSize1          5  gapopen       5   CIAGR  extended CIGAR string
#  6  name2 (query)     6  qstart        6   MRNM   Mate Reference sequence NaMe (‘=’ if same as RNAME)
#  7  start2            7  qend          7   MPOS   1-based Mate POSistion
#  8  alnSize2          8  sstart        8   ISIZE  Inferred insert SIZE
#  9  strand2           9  send          9   SEQ    query SEQuence on the same strand as the reference
#  10 seqSize2          10 evalue        10  QUAL   query QUALity (ASCII-33 gives the Phred base quality)
#  11 blocks            11 bitscore      11  OPT    variable OPTional fields in the format TAG:VTYPE:VALUE
#                                        **** OPT fields
#                                        11  AS:i:
#                                        12  XS:i:
#                                        13  XF:i:
#                                        14  XE:i:
#                                        15  NM:i:
#
#############################################################################################################

seq = {}
cov = {}
cnt = 0
length = {}

if os.path.exists(file):
    length = readFastaSeq(file)

    print("Done loading original " + str(len(length.keys())) + " sequences.")
    
with open(input_file, 'r') as input:
    for line in input:
        if line.startswith('#'):
            continue
        
        temp = line.rstrip().split('\t')
        
        if type == "guess":
            if len(temp) == 12 and temp[0].isdigit():
                type = "last"
            elif len(temp) >= 12 and temp[8].isdigit():
                type = "blast"
            elif len(temp) > 13:
                type = "sam"
            else:
                raise ValueError("ERROR: Can not recognize input format!")
            
            print("Guess input format: " + type + ".")
        
        sid, qid, qlen, qstart, qend, dist = "", "", 0, 0, 0, 1
        
        if type.lower().startswith("last"):
            if temp[6] in length:
                len = length[temp[6]]
            else:
                len = temp[10]
            
            # LAST use 0-based position
            sid = temp[1]
            qid = temp[6]
            qlen = len
            qstart = int(temp[7]) + 1
            qend = int(temp[7]) + int(temp[8])
        elif type.lower().startswith("blast"):
            if temp[0] in length:
                len = length[temp[0]]
            else:
                raise ValueError("No query sequence found.")
            
            sid = temp[1]
            qid = temp[0]
            qlen = len
            qstart = int(temp[6])
            qend = int(temp[7])
            dist = int(temp[4]) + int(temp[5])
        elif type.lower().startswith("sam"):
            nm = re.search(r'NM:i:(\d+)', line)
            if nm:
                nm = int(nm.group(1))
            length.setdefault(temp[0], len(temp[9]))
            len = length[temp[0]]
            clip5 = int(re.findall(r'^(\d+)[SH]', temp[5])[0]) if re.findall(r'^(\d+)[SH]', temp[5]) else 0
            clip3 = int(re.findall(r'(\d+)[SH]$', temp[5])[0]) if re.findall(r'(\d+)[SH]$', temp[5]) else 0
            
            if 'r' in temp[5]:
                clip5, clip3 = clip3, clip5
            
            sid = temp[2]
            qid = temp[0]
            qlen = len
            qstart = clip5 + 1
            qend = len - clip3
            dist = nm
        else:
            raise ValueError("ERROR: unknown input format.")
        
        acc = getAccFromSeqID(sid)
        taxid = acc2taxID(acc)
        length[qid] = qlen
        seq.setdefault(qid, {}).setdefault(taxid, {}).setdefault(str(qstart) + ".." + str(qend), dist)
        cov.setdefault(qid, {}).setdefault(cnt, {}).setdefault(taxid, str(qstart) + ".." + str(qend))
        
        cnt += 1

print("##SEQ\tRANK\tORGANISM\tTAX_ID\tPARENT\tLENGTH\tNUM_HIT\tTOL_HIT_LEN\tTOL_MISM\tAVG_IDT\tLINEAR_LEN\tRANK_LINEAR_LEN\tCOV\tSCALED_COV\tACC_COV_RGN\tACC_COV_LEN")

# 1  SEQ = query sequence name
# 2  RANK = rank name
# 3  ORGANISM = taxonomy name
# 4  TAX_ID = taxonomy id
# 5  PARENT = parent taxonomy name
# 6  LENGTH = query sequence name
# 7  NUM_HIT = number of hits
# 8  TOL_HIT_LEN = total bases of hits
# 9  TOL_MISM = total bases of mismatches
# 10 AVG_IDT = average hit identity
# 11 LINEAR_LEN = linear length of hits
# 12 RANK_LINEAR_LEN = total linear length of hits in the certain rank
# 13 COV = LINEAR_LEN/LENGTH
# 14 SCALED_COV = LINEAR_LEN/RANK_LINEAR_LEN

# Loop over the ranks
for rank in ["superkingdom","phylum","class","order","family","genus","species","strain"]:
    for pname in sorted(seq.keys()):
        for taxid in seq[pname].keys():
            name = taxid2rank(taxid, rank)

            # Upper taxa
            upname = taxid2rank(taxid, upper_level)
            upname = "NA" if upname is None else upname
            name = f"{upname} {rank}" if name is None else name

            if rank == "strain":
                p[rank][name] = {}
                p[rank][name]["UP_RANK"] = upname
                p[rank][name]["TAXO_ID"] = taxid
            else:
                p[rank][name] = {}
                p[rank][name]["UP_RANK"] = upname
                p[rank][name]["TAXO_ID"] = taxid2rank_taxid(taxid, rank)

            pcov = p[rank][name].get("LINEAR_LEN", "0" * length[pname])

            for region in seq[pname][taxid].keys():
                qs, qe = map(int, region.split(".."))
                end = len(pcov)
                nm = seq[pname][taxid][region]

                # Update linear length
                str_ = "0" * (qs - 1) + "1" * (qe - qs + 1) + "0" * (end - qe)
                pcov = "".join(["1" if c1 == "1" or c2 == "1" else "0" for c1, c2 in zip(pcov, str_)])
                # Total mapped
                p[rank][name].setdefault("TOL_HIT_LEN", 0)
                p[rank][name]["TOL_HIT_LEN"] += qe - qs + 1
                # Number of hits
                p[rank][name].setdefault("NUM_HIT", 0)
                p[rank][name]["NUM_HIT"] += 1
                # Distance
                p[rank][name].setdefault("TOL_MISM", 0)
                p[rank][name]["TOL_MISM"] += nm

            p[rank][name]["AVG_IDT"] = (p[rank][name]["TOL_HIT_LEN"] - p[rank][name]["TOL_MISM"]) / p[rank][name]["TOL_HIT_LEN"]
            p[rank][name]["LINEAR_LEN"] = pcov

            for name in p[rank].keys():
                p[rank][name]['LINEAR_LEN'] = p[rank][name]['LINEAR_LEN'].replace('0', '')
                sum = len(p[rank][name]['LINEAR_LEN'])

                p[rank][name]['LINEAR_LEN'] = sum
                r[rank]['TOL_LINEAR_LEN'] = r[rank].get('TOL_LINEAR_LEN', 0)
                r[rank]['TOL_LINEAR_LEN'] += sum

            # accumulated coverage
            acc_cov = '0' * length[pname]
            map = {}
            map[48] = 'unclassified'
            mid_ascii = 49

            for cnt in sorted(cov[pname].keys()):
                for taxid in cov[pname][cnt].keys():
                    name = taxid2rank(taxid, rank)
                    # upper taxa
                    upname = taxid2rank(taxid, upper_level)
                    upname = 'NA' unless upname
                    name = f'{upname} {rank}' unless name

                    if name not in map:
                        map[name] = mid_ascii
                        map[mid_ascii] = name
                        mid_ascii += 1

                    region = cov[pname][cnt][taxid]
                    acc_cov = accCov(acc_cov, chr(map[name]), region)

            csum = accCovSummary(acc_cov, map)
            for name in csum.keys():
                p[rank][name]['ACC_COV_RGN'] = ';'.join(csum[name])
                len = 0
                for rgn in csum[name]:
                    qs, qe = map(int, rgn.split('..'))
                    len += qe - qs + 1
                p[rank][name]['ACC_COV_LEN'] = len

            upper_level = rank

            for rank in ("superkingdom","phylum","class","order","family","genus","species","strain"):
                for name in sorted(p[rank], key=lambda x: p[rank][x]['ACC_COV_LEN'], reverse=True):
                    if name == "unclassified":
                        continue
                    pref = p[rank][name]
                    #sequence rank orig taxid parent #hit tol_mapped_len linear_len tol_linear_len coverage prob
                    print("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%d\t%d\t%.4f\t%.4f\t%s\t%d" % (
                        pname,
                        rank,
                        name,
                        pref['TAXO_ID'],
                        pref['UP_RANK'],
                        length[pname],
                        pref['NUM_HIT'],
                        pref['TOL_HIT_LEN'],
                        pref['TOL_MISM'],
                        pref['AVG_IDT'],
                        pref['LINEAR_LEN'],
                        r[rank]['TOL_LINEAR_LEN'],
                        pref['LINEAR_LEN']/length[pname],
                        pref['LINEAR_LEN']/r[rank]['TOL_LINEAR_LEN'],
                        pref['ACC_COV_RGN'],
                        pref['ACC_COV_LEN']
                    ))
                pref = p["superkingdom"]["unclassified"]
                print("%s\t%s\t%s\t\t\t%d\t\t\t\t\t%d\t\t%.4f\t\t%s\t%d" % (
                    pname,
                    "unclassified",
                    "unclassified",
                    length[pname],
                    pref['ACC_COV_LEN'],
                    pref['ACC_COV_LEN']/length[pname],
                    pref['ACC_COV_RGN'],
                    pref['ACC_COV_LEN']
                ))

def acc_cov(acc_cov, id, region):
    qs, qe = map(int, region.split('..'))
    while '0' in acc_cov:
        us = acc_cov.index('0') + 1
        ue = us + acc_cov[us:].index('1') - 1
        if us > qe:
            break
        if qs >= us and qe <= ue:
            # whole overlapping
            length = qe - qs + 1
            acc_cov = acc_cov[:qs-1] + id * length + acc_cov[qe:]
        elif us >= qs and qe >= us and ue >= qe:
            # cov overlapping 3" 0s
            length = ue - qs + 1
            acc_cov = acc_cov[:qs-1] + id * length + acc_cov[ue:]
        elif qs >= us and qs <= ue and qe > ue:
            # overlapping 5"
            length = qe - us + 1
            acc_cov = acc_cov[:us-1] + id * length + acc_cov[qe:]
    return acc_cov


def acc_cov_summary(acc_cov, tax_map):
    c = {}
    csum = {}
    for match in re.findall(r'(.)\1*', acc_cov):
        qs, qe = acc_cov.index(match) + 1, acc_cov.index(match) + len(match)
        tax = tax_map[ord(match)]
        prev_end = int(csum[tax][-1].split('..')[1]) if csum.get(tax) else None
        if prev_end is not None and prev_end + 1 == qs:
            csum[tax][-1] = csum[tax][-1].replace('..' + str(prev_end), '..' + str(qe))
        else:
            csum.setdefault(tax, []).append(str(qs) + '..' + str(qe))
    return csum


def read_fasta_seq(seq_file):
    hash = {}
    with open(seq_file) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                id = line[1:].split()[0]
                seq = next(fh).strip()
                length = len(seq)
                hash[id] = length
    return hash


