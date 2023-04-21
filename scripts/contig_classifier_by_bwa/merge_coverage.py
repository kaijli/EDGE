#!/usr/bin/env python

# Import required module
import sys

# Define dictionaries to hold the data
cov = {}      # Dictionary to store coverage data
length = {}   # Dictionary to store length data
old_contig = None  # Variable to store the previous contig

# Get input file and contig from command line arguments
input_file = sys.argv[1]
contig = sys.argv[2]

# Print header line
print("##SEQ\tRANK\tORGANISM\tTAX_ID\tPARENT\tLENGTH\tNUM_HIT\tTOL_HIT_LEN\tTOL_MISM\tAVG_IDT\tLINEAR_LEN\tRANK_LINEAR_LEN\tCOV\tSCALED_COV\tNUM_MERGED\tACC_COV_LEN\n")

# Open input file for reading
with open(input_file, 'r') as infile:
    # Loop through each line in the input file
    for line in infile:
        # Skip lines starting with '#' (comments)
        if line.startswith('#'):
            continue

    #SEQ    RANK    ORGANISM    TAX_ID  PARENT  LENGTH  NUM_HIT TOL_HIT_LEN TOL_MISM    AVG_IDT LINEAR_LEN  RANK_LINEAR_LEN COV SCALED_COV  ACC_COV_RGN  ACC_COV_LEN
    # 0       1        2           3      4       5        6         7         8           9        10              11       12     13           14           15
    
        # Split the line by tab ('\t') to get individual fields
        line = line.strip().split('\t')
        # Extract relevant fields from the line
        id = line[0]
        # Create nested dictionaries to store the data
        cov[contig] = {
            id: {
                line[1]: {
                    line[2]: {
                        'TAXA_ID': line[3],
                        'PARENT': line[4],
                        'LENGTH': line[5],
                        'NUM_HIT': line[6],
                        'TOL_HIT_LEN': line[7],
                        'TOL_MISM': line[8],
                        'LINEAR_LEN': line[10],
                        'ACC_COV_LEN': line[15]
                    }
                }
            }
        }
        length[contig] = {
            id: {
                'LENGTH': line[5]
            }
        }

        # Set old_contig to current contig if it's None (i.e., first iteration)
        if old_contig is None:
            old_contig = contig

        # Merge coverage data when old_contig is not the same as current contig
        if old_contig != contig:
            covergeMerger(old_contig)
            del cov[old_contig]
            del length[contig]
            old_contig = contig

# Close input file
infile.close()

# Merge coverage data for the last contig
covergeMerger(old_contig)

# 1   SEQ = query sequence name
# 2   RANK = rank name
# 3   ORGANISM = taxonomy name
# 4   TAX_ID = taxonomy id
# 5   PARENT = parent taxonomy name
# 6   LENGTH = query sequence name
# 7   NUM_HIT = number of hits
# 8   TOL_HIT_LEN = total bases of hits
# 9   TOL_MISM = total bases of mismatches
# 10  AVG_IDT = average hit identity
# 11  LINEAR_LEN = linear length of hits
# 12  RANK_LINEAR_LEN = total linear length of hits in the certain rank
# 13  COV = LINEAR_LEN/LENGTH
# 14  SCALED_COV = LINEAR_LEN/RANK_LINEAR_LEN
# 15  NUM_MERGED = number of sequences merged

def covergeMerger(contig, cov, length):
    cds_len = 0
    ctg_cov = {}
    r = {}

    for id in cov[contig]:
        cds_len += length[contig][id]['LENGTH']

        for rank in cov[contig][id]:
            for name in cov[contig][id][rank]:
                if contig not in ctg_cov:
                    ctg_cov[contig] = {}
                if rank not in ctg_cov[contig]:
                    ctg_cov[contig][rank] = {}
                if name not in ctg_cov[contig][rank]:
                    ctg_cov[contig][rank][name] = {
                        'NUM_CDS': 0,
                        'NUM_HIT': 0,
                        'TOL_HIT_LEN': 0,
                        'LINEAR_LEN': 0,
                        'TOL_MISM': 0
                    }

                cov_ctg_ref = ctg_cov[contig][rank][name]
                cov_ctg_ref['NUM_CDS'] += 1
                cov_ctg_ref['PARENT'] = cov[contig][id][rank][name]['PARENT']
                cov_ctg_ref['TAXA_ID'] = cov[contig][id][rank][name]['TAXA_ID']
                cov_ctg_ref['NUM_HIT'] += cov[contig][id][rank][name]['NUM_HIT'] or 0
                cov_ctg_ref['TOL_HIT_LEN'] += cov[contig][id][rank][name]['TOL_HIT_LEN'] or 0
                cov_ctg_ref['LINEAR_LEN'] += cov[contig][id][rank][name]['LINEAR_LEN']
                cov_ctg_ref['TOL_MISM'] += cov[contig][id][rank][name]['TOL_MISM'] or 0
                cov_ctg_ref['ACC_COV_LEN'] += cov[contig][id][rank][name]['ACC_COV_LEN']
                r[rank] = r.get(rank, 0) + cov[contig][id][rank][name]['LINEAR_LEN']

    rank_order = ["superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]
    for rank in rank_order:
        for name in sorted(ctg_cov[contig][rank], key=lambda a: ctg_cov[contig][rank][a]["ACC_COV_LEN"], reverse=True):
            cov_ctg_ref = ctg_cov[contig][rank][name]
            print("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.4f\t%d\t%d\t%.4f\t%.4f\t%d\t%d" %
                (
                  contig,                                              # 1   SEQ = query sequence name                                        
                  rank,                                                # 2   RANK = rank name                                                 
                  name,                                                # 3   ORGANISM = taxonomy name                                         
                  cov_ctg_ref["TAXA_ID"],                              # 4   TAX_ID = taxonomy id                                             
                  cov_ctg_ref["PARENT"],                               # 5   PARENT = parent taxonomy name                                    
                  cds_len,                                             # 6   LENGTH = query sequence name                                     
                  cov_ctg_ref["NUM_HIT"],                              # 7   NUM_HIT = number of hits                                         
                  cov_ctg_ref["TOL_HIT_LEN"],                          # 8   TOL_HIT_LEN = total bases of hits                                
                  cov_ctg_ref["TOL_MISM"],                             # 9   TOL_MISM = total bases of mismatches                             
                  (cov_ctg_ref["TOL_HIT_LEN"]-cov_ctg_ref["TOL_MISM"])/cov_ctg_ref["TOL_HIT_LEN"],  # 10  AVG_IDT = average hit identity                                   
                  cov_ctg_ref["LINEAR_LEN"],                           # 11  LINEAR_LEN = linear length of hits                               
                  r[rank],                                             # 12  RANK_LINEAR_LEN = total linear length of hits in the certain rank
                  cov_ctg_ref["LINEAR_LEN"]/cds_len,                    # 13  COV = LINEAR_LEN/LENGTH
                  cov_ctg_ref["LINEAR_LEN"]/r[rank],                     # 14  SCALED_COV = LINEAR_LEN/RANK_LINEAR_LEN
                  cov_ctg_ref["NUM_CDS"],                              # 15  NUM_MERGED = number of sequences merged
                  cov_ctg_ref["ACC_COV_LEN"]                           # 16  ACC_COV_LEN = number of merged ACC_COV_LEN
              )
        )

    cov_ctg_ref = ctg_cov[contig]["unclassified"]["unclassified"]
    print("%s\t%s\t%s\t\t\t%d\t\t\t\t\t%d\t\t%.4f\t\t\t%d" % (
        contig,
        "unclassified",
        "unclassified",
        cds_len,
        cov_ctg_ref['LINEAR_LEN'],
        cov_ctg_ref['LINEAR_LEN']/cds_len,
        cov_ctg_ref['ACC_COV_LEN']
    ))
