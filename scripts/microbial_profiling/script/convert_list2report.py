import sys
import getopt

opt = {}
levels_h = {}
result = {}

def usage():
    print('''USAGE: {} [--list <FILE>] [--level <RANK>] [--setting <INI_FILE>] [--top <NUM>] [--output <DIR>] [--help]

OPTIONS:

  --list <FILE>
      File containing list of files to process.

  --level <RANK>
      Default is "genus,species,strain". Any ranks separated by commas are allowable.

  --setting <INI_FILE>
      File containing settings in INI format.

  --top <NUM>
      Only display the top <NUM> results. Default is 5.

  --output <DIR>
      Directory to output results. Default is value specified in the settings file.

  --help
      Display this help message.
'''.format(sys.argv[0]))
    sys.exit()

# Parse command line options
try:
    opts, args = getopt.getopt(sys.argv[1:], "l:s:t:o:h", ["list=", "level=", "setting=", "top=", "output=", "help"])
except getopt.GetoptError as err:
    print(str(err))
    usage()

for opt, arg in opts:
    if opt in ("-l", "--level"):
        levels_h = {lvl:1 for lvl in arg.split(",")}
    elif opt in ("-s", "--setting"):
        opt['setting'] = arg
    elif opt in ("-t", "--top"):
        opt['top'] = int(arg)
    elif opt in ("-o", "--output"):
        opt['output'] = arg
    elif opt in ("-h", "--help"):
        usage()

opt['level'] = "genus,species,strain" if 'level' not in opt else opt['level']
opt['top'] = 5 if 'top' not in opt else opt['top']

# Read settings
ini_file = opt['setting']
tools = restore_settings(ini_file)

file_info = restore_filelist(opt['list'])

# Default settings
p_outdir = tools['system']['OUTDIR']
p_outdir = opt['output'] if 'output' in opt else p_outdir

p_repdir = p_outdir + "/" + tools['system']['REPDIR']
count = 1

# Submit tools
for idx in sorted(file_info.keys()):
    fnb = file_info[idx]['PREFIX']

    for tool in sorted(tools.keys(), key=lambda x: tools[x]['ORDER']):
        if tool == 'system':
            continue
        file = "{}/{}/{}-{}/{}-{}.list.txt".format(p_repdir, count, fnb, tool, fnb, tool)
        print("Processing {}...".format(file), , file=sys.stderr)

        with open(file, 'r') as f:
            for line in f:
                fields = line.strip().split("\t")
                lvl = fields[0]
                if lvl in levels_h and fields[2]:
                    result[fnb][tool][lvl][fields[1]] = fields[2]

    count += 1

print("Creating TOP{} list...".format(opt['top']), file=sys.stderr)

# Print header
print("DATASET\tTOOL\tLEVEL", end='')
for i in range(1, opt['top'] + 1):
    print("\tTOP{}".format(i), end='')
print()

# Print content
for idx in sorted(file_info.keys()):
    fnb = file_info[idx]['PREFIX']
    for tool in sorted(tools.keys(), key=lambda x: tools[x]['ORDER']):
        if tool == 'system':
            continue
        for lvl in levels_a:
            print("{}\t{}\t{}".format(fnb, tool, lvl), end='')
            
            # No results at certain level
            if result.get(fnb, {}).get(tool, {}).get(lvl) is None:
                print("\tN/A" * opt['top'])
                continue
            
            # Print top # results
            count = 0
            for taxa in sorted(result[fnb][tool][lvl].keys(), key=lambda x: (result[fnb][tool][lvl][x], x), reverse=True):
                print("\t{}".format(taxa), end='')
                count += 1
                if count == opt['top']:
                    break
            
            # Fill "N/A" if the results do not fulfill top number
            print("\tN/A" * (opt['top'] - count))
            print()

def restore_settings(file):
    set_data = {}
    with open(file, 'r') as f:
        count = 0
        section = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#') or line.startswith(';'):
                continue
            if not line.strip():
                continue
            if line.startswith('[') and line.endswith(']'):
                section = line[1:-1]
                set_data[section] = {}
                set_data[section]['ORDER'] = count
                count += 1
                continue
            key, val = line.split('=', 1)
            set_data[section][key] = val
    return set_data

def restore_filelist(file):
    file_info = {}
    with open(file, 'r') as f:
        count = 1
        filelist_header = None
        for line in f:
            line = line.strip()
            if line.startswith('--'):
                continue
            if not line:
                continue
            if line.startswith('#'):
                continue
            if line.startswith('PREFIX'):
                filelist_header = line.split('\t')
                continue
            fields = line.split('\t')
            for i in range(len(filelist_header)):
                file_info[count][filelist_header[i]] = fields[i]
            count += 1
    return file_info
