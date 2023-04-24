import sys
import getopt
import random

def usage():
    print("python3 random_seq_extractor.py -o <output> -i <fasta/q> -n #")
    print("     -i    input fasta/q file")
    print("     -o    output fasta/q file")
    print("     -n    number of sequences to extract [int]")
    print("           it can be a fraction of total sequences [0.0 to 1.0]")
    print("     -h    print usage")
    sys.exit()

def seq_count(file):
    count = 0
    mode = "fasta"
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('@'):
                mode = "fastq"
            count += 1
    if mode == "fastq":
        count = count // 4
    return count

def main():
    output = ''
    file = ''
    output_num = 0
    try:
        opts, args = getopt.getopt(sys.argv[1:], "o:i:n:h")
    except getopt.GetoptError as err:
        print(str(err))
        usage()
    for opt, arg in opts:
        if opt == '-o':
            output = arg
        elif opt == '-i':
            file = arg
        elif opt == '-n':
            output_num = float(arg)
        elif opt == '-h':
            usage()
    if not (output and file and output_num):
        usage()

    total_seq_num = seq_count(file)
    if output_num >= total_seq_num:
        print("The extract number({}) is bigger/equal than sequence total number({}).".format(output_num, total_seq_num))
        sys.exit(1)

    if output_num > 0 and output_num < 1:
        output_num = int(output_num * total_seq_num)

    output_num = int(output_num)
    print("Extracting random {} sequences out of total {}".format(output_num, total_seq_num))

    get = {}
    while True:
        get[random.randint(1, total_seq_num)] = 1
        num_of_random = len(get)
        if num_of_random == output_num:
            break
    print("Done choosing {} sequences.".format(output_num))

    name = ''
    fastq = False
    count = 0

    with open(output, 'w') as out, open(file, 'r') as inp:
        for line in inp:
            line = line.strip()
            if line.startswith('@'):
                count += 1
                name = line
                seq = inp.readline().strip()
                while not seq.startswith('+'):
                    seq += inp.readline().strip()
                q_name_pos = seq.index('+')
                qual_id = seq[q_name_pos:]
                seq = seq[:q_name_pos]
                seq_len = len(seq)
                qual_seq = inp.readline().strip()
                qual_seq_len = len(qual_seq)
                while qual_seq_len < seq_len:
                    if qual_seq_len == seq_len:
                        break
                    qual_seq += inp.readline().strip()
                    qual_seq_len = len(qual_seq)
                print_out = "{}\n{}\n{}\n{}".format(name, seq, qual_id, qual_seq)
                if count in get:
                    out.write(print_out + '\n')
            elif line.startswith('>'):
                count += 1
                if count in get:
                    out.write(line + '\n')
            else:
                if count in get:
                    out.write(line + '\n')

if __name__ == '__main__':
    main()
