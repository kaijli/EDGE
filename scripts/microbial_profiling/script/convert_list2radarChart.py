#!/usr/bin/python3
import os
import re
import sys
import getopt
from pathlib import Path
from htmltemplate import Template

template = None
opt = {}
res, args = getopt.getopt(sys.argv[1:], 'l:o:p:f:t:h', ['level=', 'outdir=', 'outprefix=', 'filter_taxa=', 'title=', 'template=', 'top=', 'help'])

def usage():
    print("Usage: script.py [options] <listfile(s)>")
    print("Options:")
    print("  -l, --level=<s>         : Level (default: genus)")
    print("  -o, --outdir=<s>        : Output directory (default: .)")
    print("  -p, --outprefix=<s>     : Output file prefix (default: \"\")")
    print("  -f, --filter_taxa=<s>   : Filter taxa (default: \"\")")
    print("  -t, --template=<s>      : Template file (default: \"convert_list2radarChart.html.tmpl\")")
    print("  --top=<i>               : Top N (default: 5)")
    print("  -h, --help              : Show this help message")
    sys.exit(1)

for option, value in res:
    if option in ['-l', '--level']:
        opt['level'] = value
    elif option in ['-o', '--outdir']:
        opt['outdir'] = value
    elif option in ['-p', '--outprefix']:
        opt['outprefix'] = value
    elif option in ['-f', '--filter_taxa']:
        opt['filter_taxa'] = value
    elif option in ['-t', '--template']:
        opt['template'] = value
    elif option == '--top':
        opt['top'] = int(value)
    elif option in ['-h', '--help']:
        usage()

if not args:
    usage()

opt['template'] = opt.get('template', os.path.join(os.path.dirname(__file__), '../../../lib/convert_list2radarChart.html.tmpl'))
opt['outdir'] = opt.get('outdir', '.')
opt['outprefix'] = opt.get('outprefix', '')
opt['top'] = opt.get('top', 5)
opt['filter_taxa'] = opt.get('filter_taxa', '')
opt['level'] = opt.get('level', 'genus')

if Path(opt['template']).is_file():
    template = Template(filename=opt['template'])
else:
    print(f"ERROR: Can't find template: {opt['template']}.")
    sys.exit(1)

filter_taxa = set(opt['filter_taxa'].split(','))
listfile = args
matrix = None
count = None
order = None
data = None
tools = None
tool_order = {
	'blastn': 10,
	'bwa': 20,
	'bwa-target': 21,
	'metaphlan': 22,
	'metaphlan2': 23,
	'motus': 24,
	'kraken_mini': 25,
	'kraken': 26,
	'kraken2': 27,
	'centrifuge': 28,
	'gottcha-genDB-b': 30,
	'gottcha-speDB-b': 40,
	'gottcha2-speDB-b': 41,
	'gottcha-strDB-b': 50,
	'gottcha-genDB-v': 60,
	'gottcha-speDB-v': 70,
	'gottcha2-speDB-v': 71,
	'gottcha-strDB-v': 80,
	'sequedex-opg': 90,
	'sequedex-tol': 100,
	'metacv': 120,
	'metaphyler-bn': 130,
	'metaphyler-bx': 140,
	'metaphyler-srv': 150,
	'pangia': 160,
	'diamond': 170
}

############################
# parsing list files
############################

for file in listfile:
    dataset, tool = re.findall(r'\d+_([^\/]+)\/([^\/]+)\/([^\/]+)\.[^\.]+$', file)[0]

    if 'species' in opt['level']:
        if re.search(r'gottcha-genDB', tool):
            continue
    elif 'strain' in opt['level']:
        if re.search(r'gottcha-genDB', tool) or re.search(r'gottcha-speDB', tool):
            continue

    tools[tool] = 1
    filtered_taxa = 0

    with open(file, 'r') as list_file:
        for line in list_file:
            fields = line.strip().split('\t')
            if fields[0] != opt['level']:
                continue
            if fields[1] in filter_taxa:
                continue

            matrix.setdefault(dataset, {}).setdefault(tool, {})[fields[1]] = fields[2]

            # counting total reads/abundance
            count.setdefault(dataset, {}).setdefault(tool, 0)
            count[dataset][tool] += fields[2]

#############################
# Create orders
#############################

top_taxa = {}

for dataset in matrix.keys():
    for tool in matrix[dataset].keys():
        cnt = 0
        top_sum = 0
        for taxa in sorted(matrix[dataset][tool], key=lambda x: matrix[dataset][tool][x], reverse=True):
            if cnt == opt['top']:
                break
            top_taxa[dataset][taxa] = 1
            top_sum += matrix[dataset][tool][taxa]
            cnt += 1
        matrix[dataset][tool]['others'] = count[dataset][tool] - top_sum

#############################
# Generating matrix
#############################

legend_text = ""
data_text = ""

for dataset in top_taxa.keys():
    legend = sorted(tools.keys(), key=lambda x: tool_order[x])
    data_tool = []

    print("Generate chart: ", dataset)

    for tool in sorted(tools.keys(), key=lambda x: tool_order[x]):
        data = []
        top_taxas = sorted(top_taxa[dataset].keys(), key=lambda x: (matrix[dataset]['bwa'][x], matrix[dataset]['gottcha-genDB-b'][x]))
        for taxa in ["others"] + top_taxas:
            val = matrix[dataset][tool].get(taxa, 0)
            tol = count[dataset][tool]
            tol = 1 if val == 0 and tol == 0 else tol
            pct = round(val/tol, 4)
            data.append("{{'axis':'{}', 'value': {}}}".format(taxa, pct))
            print("\t{}\t{}\t{}\t{}".format(tool, taxa, val, pct))
        temp = ",\n".join(data)
        data_tool.append("[ {} ]".format(temp))

    legend_text = "','".join(legend)
    legend_text = "'{}'".format(legend_text)

    data_text = ",\n".join(data_tool)

    template.param('TITLE', opt['title'] if 'title' in opt else dataset)
    template.param('LEGEND', legend_text)
    template.param('DATA', data_text)

    outprefix = opt['outprefix'] + "_" if 'outprefix' in opt else ""
    with open("{}/{}{}.{}.html".format(opt['outdir'], outprefix, dataset, opt['level']), "w") as OUTPUT:
        OUTPUT.write(template.output())

import sys

def usage():
    print('''USAGE: {} <DATASET>/<TOOL>/<DATASET>.out.list [OPTIONS]

EXAMPLE: 
    {} 1_CDC_10_clean/bwa/*.out.list 1_CDC_10_clean/metaphlan/*.out.list

    {} \\
        --level genus \\
        --top 5 \\
        --template convert_list2radarChart.html.tmpl \\
        1_454_EVEN/*/*.out.list

OPTIONS:

  --level <RANK>

      Default is "species". Any ranks from "kingdom" to "strain" are allowable.

  --top (INT)

      Only display a certain number of classifications.

  --exclude <TAXA1(,TAXA2,...)>
  --outdir <DIR>
  --outprefix <STRING>
  --template <HTML TEMPLATE>

  --help
'''.format(sys.argv[0], sys.argv[0], sys.argv[0]))
    sys.exit()

