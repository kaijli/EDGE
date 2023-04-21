#!/usr/bin/env python
# contig_classifier.py
# ver 0.1
# 2013/09/25
#
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.

# Change log
# 2015 Feb
# - added accumulated coverage
# - added LCA

import os
import sys
import subprocess
import time
import getopt

def usage():
    """Print usage information"""
    print("""Usage: contig_classifier.py -i <input_file> [-p <prefix>] [-t <threads>] [-d <database>] [--pacbio] [--debug] [--help]
Options:
    -i, --input <input_file>      Input file in SAM format
    -p, --prefix <prefix>        Prefix for output files (default: 'output')
    -t, --threads <threads>      Number of threads to use (default: 2)
    -d, --db <database>          Database file for classification (default: '/opt/apps/edge/database/bwa_index/NCBI-Bacteria-Virus.fna')
        --pacbio                  Use PacBio mode for BWA (default: off)
        --debug                    Print debug information
    -h, --help                    Print this help message
""")

def execute_command(cmd):
    """Execute a command in the shell"""
    subprocess.run(cmd, shell=True, check=True)

def time_interval(start_time):
    """Calculate time interval since start time"""
    end_time = time.time()
    interval = end_time - start_time
    return time.strftime("%H:%M:%S", time.gmtime(interval))

def mem_usage():
    """Get memory usage"""
    process = subprocess.Popen("ps -p {} -o rss | tail -1".format(os.getpid()), stdout=subprocess.PIPE, shell=True)
    output, _ = process.communicate()
    return int(output.strip()) / 1024

def count_result(file):
    """Count the number of contigs and bases in the result file"""
    contigs_count = 0
    contigs_bases = 0
    classified_contigs_count = 0
    classified_contigs_bases = 0
    unclassified_contigs_count = 0
    unclassified_contigs_bases = 0
    with open(file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) == 14:
                contigs_count += 1
                contigs_bases += int(fields[5])
                if fields[3] != 'unclassified':
                    classified_contigs_count += 1
                    classified_contigs_bases += int(fields[5])
                else:
                    unclassified_contigs_count += 1
                    unclassified_contigs_bases += int(fields[5])
    return contigs_count, contigs_bases, classified_contigs_count, classified_contigs_bases, unclassified_contigs_count, unclassified_contigs_bases

def main(argv):
    input_file = ''
    pacbio = False
    threads = 2
    prefix = 'output'
    db = '/opt/apps/edge/database/bwa_index/NCBI-Bacteria-Virus.fna'
    debug = False
    help_flag = False

    try:
        opts, args = getopt.getopt(argv, "hi:t:p:d:", ["input=", "pacbio", "threads=", "prefix=", "db=", "debug", "help"])
    except getopt.GetoptError as e:
        print(str(e))
        usage()

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            help_flag = True
        elif opt in ("-i", "--input"):
            input_file = arg
        elif opt == "--pacbio":
            pacbio = True
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-p", "--prefix"):
            prefix = arg
        elif opt in ("-d", "--db"):
            db = arg
        elif opt == "--debug":
            debug = True

    if help_flag or not input_file:
        usage()

    filename = os.path.basename(input_file)
    basename, ext = os.path.splitext(filename)

    # Initialize options
    threads = threads or 2
    prefix = prefix or 'output'

    start_time = time.time()
    period = time_interval(start_time)
    mem = mem_usage()
    period_str = f"[{period}]"

    mem = mem_usage()
    print(f"{period_str} Running BWA")

    if opt["pacbio"]:
        execute_command(f"bwa mem -B5 -Q2 -E1 -a -M -t {threads} {db} temp$$/{filename} > temp$$/{filename}.sam 2>>{prefix}.log")
    else:
        execute_command(f"bwa mem -a -M -t {threads} {db} temp$$/{filename} > temp$$/{filename}.sam 2>>{prefix}.log")

    period = time_interval(time)
    print(f"{period_str} Filtering unmapped contigs")
    execute_command(f"samtools view -F4 temp$$/{filename}.sam > temp$$/{filename}.mapped.sam 2>>{prefix}.log")

    period = time_interval(time)
    print(f"{period_str} Splitting SAM file")
    execute_command(f"split_sam_by_lines.pl --line 20000 --input temp$$/{filename}.mapped.sam")

    period = time_interval(time)
    print(f"{period_str} Classifying contigs")
    execute_command(f"cd temp$$; parallel --results {prefix}_para_log -j {threads} 'seq_coverage.pl --input {{}}' ::: *.part* > ../{prefix}.ctg_class.csv 2>>../{prefix}.log")

    period = time_interval(time)
    print(f"{period_str} Reporting unclassified contigs")
    execute_command(f"samtools view -f4 -S temp$$/{filename}.sam 2>>{prefix}.log | awk -F\\\\t '{{print $1\"\\tunclassified\\tunclassified\\t\\t\\t\",length($10),\"\\t\\t\\t\\t\\t0\\t\\t0\\t\\t\\t\",length($10)}}' >> {prefix}.ctg_class.csv")
    execute_command(f"samtools view -f4 -S temp$$/{filename}.sam 2>>{prefix}.log | awk '{{print \">\"$1\"\\n\"$10}}' > {prefix}.unclassified.fasta")

    period = time_interval(time)
    print(f"{period_str} Merging classification")
    execute_command(f"merge_coverage.pl {prefix}.ctg_class.csv {prefix} > {prefix}.assembly_class.csv 2>>{prefix}.log")

    period = time_interval(time)
    print(f"{period_str} Reporting classification (BEST hit)")
    execute_command(f"class_top_hit_summary.pl < {prefix}.assembly_class.csv > {prefix}.assembly_class.top.csv 2>>{prefix}.log")
    execute_command(f"class_top_hit_summary.pl <

    period = time_interval(time_start)
    print(f"[{period}] Reporting classification (LCA)")
    execute_command(f"report_LCA.pl < {PREFIX}.ctg_class.csv > {PREFIX}.ctg_class.LCA.csv 2>>{PREFIX}.log")
    execute_command(f"(head -n 1 {PREFIX}.ctg_class.LCA.csv && tail -n +2 {PREFIX}.ctg_class.LCA.csv | sort -t '\t' -k6nr) > {PREFIX}.ctg_class.LCA.csv.sort")
    execute_command(f"mv {PREFIX}.ctg_class.LCA.csv.sort {PREFIX}.ctg_class.LCA.csv")

    tol_contigs_count, tol_contigs_bases, classified_contigs_count, classified_contigs_bases, unclassified_contigs_count, unclassified_contigs_bases = count_result(f"{PREFIX}.ctg_class.top.csv")
    period = time_interval(time_start)
    print(f"[{period}] Total Contigs: {tol_contigs_count} ({tol_contigs_bases} bp); Classified Contigs: {classified_contigs_count} ({classified_contigs_bases} bp); Unclassified Contigs: {unclassified_contigs_count} ({unclassified_contigs_bases} bp);")

    if not opt['debug']:
        period = time_interval(time_start)
        print(f"[{period}] Cleaning temporary files")
        execute_command(f"rm -rf temp$$/")

period = time_interval(time_start)
print(f"[{period}] Finished. Please find the result for contig in {PREFIX}.ctg_class.csv and {PREFIX}.assembly_class.csv")
