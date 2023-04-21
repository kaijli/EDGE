#!/usr/bin/env python
# Original Perl script with comments retained

acc = {}
name = ""

rank = {
    "unclassified": 0,
    "superkingdom": 1,
    "phylum": 2,
    "class": 3,
    "order": 4,
    "family": 5,
    "genus": 6,
    "species": 7,
    "strain": 8
}

print("##SEQ\tRANK\tORGANISM\tTAX_ID\tPARENT\tLENGTH\tNUM_HIT\tTOL_HIT_LEN\tTOL_MISM\tAVG_IDT\tLINEAR_LEN\tRANK_LINEAR_LEN\tCOV\tSCALED_COV\tNUM_MERGED\tACC_COV_LEN")

for line in iter(input, ''):
    if line.startswith('#'):
        continue
    line = line.rstrip()
    if "\t0$" in line:
        continue
    temp = line.split('\t')
    # SEQ    RANK    ORGANISM    TAX_ID  PARENT  LENGTH  NUM_HIT TOL_HIT_LEN TOL_MISM    AVG_IDT LINEAR_LEN  RANK_LINEAR_LEN COV SCALED_COV  ACC_COV_RGN  ACC_COV_LEN
    # 0       1        2           3      4       5        6         7         8           9        10              11       12     13           14           15

    if temp[0] != name:
        for r in sorted(rank.keys(), key=lambda x: rank[x], reverse=True):
            if r not in acc:
                continue
            if len(acc[r]) == 1:
                print(acc[r][0])
                break
        acc = {}
        name = temp[0]

    acc.setdefault(temp[1], []).append(line)
