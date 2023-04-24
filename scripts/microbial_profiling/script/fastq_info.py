#!/usr/bin/env python3
import sys
import getopt

inFastq = ''
maxReads = 0

try:
    opts, args = getopt.getopt(sys.argv[1:], "i:m:", ["i=", "m="])
except getopt.GetoptError as err:
    print(str(err))
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-i", "--i"):
        inFastq = arg
    elif opt in ("-m", "--m"):
        maxReads = int(arg)

if not inFastq:
    print("ERROR: Input FASTQ file not specified.")
    sys.exit(2)

bases = 0
reads = 0
header = ""
min_score = 164
max_score = 0
tol_avg_score = 0
tol_read_num = 0

if maxReads > 0:
    with open(inFastq) as f:
        tol_read_num = sum(1 for line in f) // 4

with open(inFastq) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith('@'):
            header = line
        elif header:
            seq = line
            bases += len(seq)
            reads += 1
            if reads % 5000 == 0:
                print(f"processing {reads} reads in {inFastq}\r", end='', file=sys.stderr)

            tmp = next(f)  # +
            qual_seq = next(f).strip()

            # score
            qual = [ord(ascii) for ascii in qual_seq]
            score = sum(qual)
            min_score = min(min_score, *qual)
            max_score = max(max_score, *qual)

            tol_avg_score += score / len(seq)

            if maxReads and reads >= maxReads:
                break

tol_read_num = reads if not maxReads else tol_read_num

platform = guessPlatform(header)
format, offset = guessOffset(min_score, max_score)


try:
    opts, args = getopt.getopt(sys.argv[1:], "i:m:", ["inFastq=", "maxReads="])
except getopt.GetoptError as err:
    print(str(err))
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-i", "--inFastq"):
        inFastq = arg
    elif opt in ("-m", "--maxReads"):
        maxReads = int(arg)

if not inFastq:
    print("Usage: python script.py -i <fastq file> [-m (PROCESS_NUMBER_OF_READS)]")
    sys.exit(2)



if maxReads > 0:
    with open(inFastq, 'r') as f:
        tol_read_num = sum(1 for _ in f) // 4

with open(inFastq, 'r') as f:
    for line in f:
        line = line.strip()
        if reads % 5000 == 0:
            print("processing {} reads in {}\r".format(reads, inFastq), end='', file=sys.stderr)

        if header == '':
            header = line
        elif header != '' and len(header) > 0:
            seq = line
            bases += len(seq)
            reads += 1

            tmp = f.readline().strip()  # +

            qual_seq = f.readline().strip()
            qual = list(qual_seq)
            score = sum(ord(ascii) for ascii in qual)
            min = min(ord(ascii), min)
            max = max(ord(ascii), max)

            tol_avg_score += score / len(seq)

            if maxReads > 0 and reads >= maxReads:
                break

    tol_read_num = reads if maxReads == 0 else tol_read_num

platform = guessPlatform(header)
format, offset = guessOffset(min, max)

print("{:14s}{:14s}{:16s}{:15s}{:14s}{:11s}{:8s}{:18s}".format("TOL_READS", "PROCESSED", "PLATFORM", "TOL_BASES", "AVG_LENGTH", "AVG_SCORE", "OFFSET", "FASTQ_FMT"))
print("----------------------------------------------------------------------------------------------------------------")
print("{:14s}{:14s}{:16s}{:15s}{:14.2f}{:11.2f}{:8d}{:18s}".format(str(tol_read_num), str(reads), platform, str(bases), (bases/reads), (tol_avg_score/reads-offset), offset, format))


def guessPlatform(header):
    header = header[1:]  # remove leading '@'
    guess_ilu = header.split(':')
    guess_pac = header.split('_')
    guess = []

    if "MiSeq" in header or "HiSeq" in header:
        guess.append("illumina")
    elif len(guess_ilu) >= 5:
        guess.append("Illumina")
    elif len(guess_pac) >= 5 and header.startswith('m.*_s\d+_p\d+'):
        guess.append("pacbio")
    elif len(header) == 14:
        guess.append("ionTorrent")
    else:
        guess.append("unknown")

    return ' '.join(guess)


def guessOffset(min, max):
    if min < 59:
        if max == 74:
            return "Illumina 1.8+", 33
        else:
            return "Sanger", 33
    elif max > 74:
        if min >= 67:
            return "Illumina 1.5+", 64
        elif min >= 64:
            return "Illumina 1.3+", 64
        else:
            return "Solexa", 64
    else:
        return "Unknown", "NA"

def usage():
    print("Usage: $0 -i <fastq file> [-m (PROCESS_NUMBER_OF__READS)]")
    exit()

# Note: 'exit()' function in Perl is equivalent to 'exit()' function in Python
