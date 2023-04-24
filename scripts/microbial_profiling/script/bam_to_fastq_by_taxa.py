#!/usr/bin/env python
# require samtools 
import os
import sys
import getopt
import time
from gi2lineage import *

workingDir = os.getcwd()
in_offset = 33
out_offset = 33
opt_mapped = False
opt_unmapped = False
opt_passfilter = False
prefix = "Reads"
mapped_ref_id = "" # this is the id in SAM 'RNAME' field (3rd column)
sam_format = False
rank = "genus"
mapped_ref_name = "" # this is the name in a taxonomy rank
preload = False
single_end = False
zip = False
fastq_file = ""
reads_id = {}

def Usage():
    print("Usage: python script.py [options] file")
    # Add usage information here

try:
    opts, args = getopt.getopt(sys.argv[1:], "S", ["in_offset=", "out_offset=", "mapped", "id=", "rank=", "name=", "fastq=", "preload", "unmapped", "pf", "se", "zip", "prefix=", "help"])
except getopt.GetoptError as err:
    print(str(err))
    Usage()
    sys.exit(2)

time = int(time.time())
if not args:
    Usage()

for opt, arg in opts:
    if opt in ("-S", "--sam"):
        sam_format = True
    elif opt == "--in_offset":
        in_offset = int(arg)
    elif opt == "--out_offset":
        out_offset = int(arg)
    elif opt == "--mapped":
        opt_mapped = True
    elif opt == "--id":
        mapped_ref_id = arg
    elif opt == "--rank":
        rank = arg
    elif opt == "--name":
        mapped_ref_name = arg
    elif opt == "--fastq":
        fastq_file = arg
    elif opt == "--preload":
        preload = True
    elif opt == "--unmapped":
        opt_unmapped = True
    elif opt == "--pf":
        opt_passfilter = True
    elif opt == "--se":
        single_end = True
    elif opt == "--zip":
        zip = True
    elif opt == "--prefix":
        prefix = arg
    elif opt in ("-h", "--help"):
        Usage()
        sys.exit()

if in_offset and not out_offset:
    sys.exit("Please provide output offset for conversion.")
if not in_offset and out_offset:
    sys.exit("Please provide input offset for conversion.")
if opt_mapped and opt_unmapped:
    sys.exit("Please provide either --mapped or --unmapped.")

file = args[0]
samtools_flag = ""
if opt_mapped or mapped_ref_name or mapped_ref_id:
    samtools_flag = "-F 4"
if opt_unmapped:
    samtools_flag = "-f 4"

with open(prefix + ".1.fastq", "w") as pair1_fh, open(prefix + ".2.fastq", "w") as pair2_fh, open(prefix + ".se.fastq", "w") as se_fh:
    sortName_sam_file = "unmapped" + str(time) + ".sam"
    if mapped_ref_id:
        mapped_ref_id = "\"" + mapped_ref_id + "\""
    if sam_format:
        samtools_flag += " -S"

    # Convert preload variable to string 'preload' if True
    if preload:
        preload = 'preload'

    if mapped_ref_name:
        print_timeInterval(time, "Load Tanonomy info ... ")
        loadTaxonomy(preload)

    if not single_end:
        print_timeInterval(time, "Sorting the SAM file ... ")
        os.system(f"samtools view {samtools_flag} {file} {mapped_ref_id} | sort -T {workingDir} -k 1,1 > {sortName_sam_file}")
        print_timeInterval(time, "Parsing the SAM file ...")
        with open(sortName_sam_file, "r") as IN:
            # continue with the rest of the code
    else:  # skip sorting to speed up process
        print_timeInterval(time, "Parsing the SAM file ...")
        IN = os.popen(f"samtools view {samtools_flag} {file} {mapped_ref_id}")
        # continue with the rest of the code

    p1_count = 0
    p2_count = 0
    se_count = 0

    # Open input and output files
    with open('IN', 'r') as in_file:
        for line in in_file:
            line = line.rstrip('\n')
            if line.startswith('@SQ'):
                continue
            array = line.split('\t')
            if opt_passfilter and (int(array[1]) & 512):
                continue
            if mapped_ref_name:
                acc = get_acc_from_seq_id(array[2])
                acc_to_name = acc2rank(acc, rank.lower())
                if acc_to_name.lower() != mapped_ref_name.lower():
                    continue
            if fastq_file:
                reads_id[array[0]] = 1
                tmp = array[0].rsplit('_', 1)[0]
                reads_id[tmp] = 1
                continue
            if in_offset != out_offset:
                array[10] = quality_conversion(array[10], in_offset, out_offset)
            if int(array[1]) & 1:
                if (opt_mapped and int(array[1]) & 8) or \
                        (opt_unmapped and not (int(array[1]) & 4) and not (int(array[1]) & 8)):
                    se_count += 1
                    if int(array[1]) & 16:
                        array[9] = reverse_complement(array[9])
                        array[10] = array[10][::-1]
                    se_fh.write('@{}\n'.format(array[0]))
                    se_fh.write('{}\n'.format(array[9]))
                    se_fh.write('+\n')
                    se_fh.write('{}\n'.format(array[10]))
                else:
                    if int(mapped_ref_id) and array[6] != '=':
                        se_count += 1
                        se_fh.write('@{}\n'.format(array[0]))
                        se_fh.write('{}\n'.format(array[9]))
                        se_fh.write('+\n')
                        se_fh.write('{}\n'.format(array[10]))
                        continue
                    if int(array[1]) & 64:
                        p1_count += 1
                        pair1_fh.write('@{}/1\n'.format(array[0]))
                        pair1_fh.write('{}\n'.format(array[9]))
                        pair1_fh.write('+\n')
                        pair1_fh.write('{}\n'.format(array[10]))
                    elif int(array[1]) & 128:
                        p2_count += 1
                        pair2_fh.write('@{}/2\n'.format(array[0]))
                        pair2_fh.write('{}\n'.format(array[9]))
                        pair2_fh.write('+\n')
                        pair2_fh.write('{}\n'.format(array[10]))
            else:
                se_count += 1
                if int(array[1]) & 16:
                    array[9] = reverse_complement(array[9])
                    array[10] = array[10][::-1]
                se_fh.write('@{}\n'.format(array[0]))
                se_fh.write('{}\n'.format(array[9]))
                se_fh.write('+\n')
                se_fh.write('{}\n'.format(array[10]))

    in_file.close()

    if fastq_file:
        extract_from_original_fastq(fastq_file, reads_id)

pair1_fh.close()
pair2_fh.close()
se_fh.close()

if rank:
    print(f"{rank}: ")
if mapped_ref_name:
    print(mapped_ref_name)
print(f"Paired 1: {p1_count}; ", end="")
print(f"Paired 2: {p2_count}; ", end="")
print(f"Single End: {se_count}; ")

if not os.path.getsize(f"{prefix}.1.fastq"):
    os.unlink(f"{prefix}.1.fastq")
if not os.path.getsize(f"{prefix}.2.fastq"):
    os.unlink(f"{prefix}.2.fastq")
if not os.path.getsize(f"{prefix}.se.fastq"):
    os.unlink(f"{prefix}.se.fastq")
if not single_end:
    os.unlink(sortName_sam_file)

if zip:
    print_timeInterval(time, "Compressed fastq files ... \n")
    file_name, file_path, file_suffix = fileparse(prefix, r"\.[^.]*")
    os.chdir(file_path)
    shutil.make_archive(file_name, "zip", ".", file_name + "*fastq")
    for file in glob.glob(file_name + "*fastq"):
        os.unlink(file)

print_timeInterval(time, "All Done.\n")

def extract_from_original_fastq(fastq, reads_id):
    with open(fastq, "r") as fh:
        for line in fh:
            s_id = line.strip()
            seq = next(fh).strip()
            q_id = next(fh).strip()
            q_seq = next(fh).strip()
            if s_id.startswith("@"):
                read_id = s_id[1:]
                if read_id in reads_id:
                    se_fh.write(f"{s_id}\n{seq}\n{q_id}\n{q_seq}\n")
                    global se_count
                    se_count += 1


def quality_conversion(seq, in_offset, out_offset):
    seq_len = len(seq)
    converted_seq = ""
    for i in range(seq_len):
        q = seq[i]
        cov_q = chr(ord(q) - in_offset + out_offset)
        converted_seq += cov_q
    return converted_seq

def reverse_complement(dna):
    reverse_comp_seq = dna[::-1]
    trans = str.maketrans("atgcrywsmkATGCRYWSMK", "tacgyrswkmTACGYRSWKM")
    reverse_comp_seq = reverse_comp_seq.translate(trans)
    return reverse_comp_seq


def usage():
    print("""Usage: <bam_file>
# require samtools in PATH.

Options:
-sam|S          Input is sam format.
-in_offset      Offset number for ASCII encoding from input
-out_offset     Offset number for ASCII encoding for output
                Type of ASCII encoding: 33 (standard) or 64 (illumina 1.3+)
                Use this option to convert the encoding offset.
                default: no conversion.
-mapped         For retrieving mapped reads from alignment bam file
                -id   retrieving reads mapped to this reference ID. (e.g. contig_00001)
                  Or   retrieving reads mapped to a specified taxonomy name
                -rank  phylum, class, order, family, genus, species or strain (default: genus)
                -name  the name in a taxonomy rank (e.g. Pseudomonas)

-unmapped       For retrieving unmapped reads from alignment bam file
-pf             Output passFilter reads only
-se             The sam file is from single end reads, will bypass sorting to speed up.
-zip            Compress output files with tar gz.
-fastq          Original mapping fastq file. Will use this to extract reads for single end reads only
-pefix          Output file prefix (default: Reads)
-help           Show this usage

""")
    #-type           pe or se.  pe=paired end, se=single end. The pe type will output reads' name with /1 or /2.
    exit()


def print_time_interval(now, msg):
    now = time.time() - now
    string = time.strftime("%H:%M:%S", time.gmtime(now))
    print(f"[{string}] {msg}\n")

