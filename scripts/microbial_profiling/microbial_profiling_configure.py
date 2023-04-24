#!/usr/bin/env python
import sys
import getopt
import json
from os.path import dirname
sys.path.append(dirname(dirname(__file__)) + '/lib')
from htmltemplate import Template

# setting up default values
opt = {
    'template': sys.argv[1] if len(sys.argv) > 1 else '',
    'tools': sys.argv[2] if len(sys.argv) > 2 else '',
    'bwaScoreCut': 30,
    'splitrim-minq': 20,
    'bwa-db': dirname(dirname(__file__)) + '/database/bwa_index/NCBI-Bacteria-Virus.fna',
    'metaphlan-db': dirname(dirname(__file__)) + '/database/metaphlan2',
    'kraken-db': dirname(dirname(__file__)) + '/database/kraken2/taxo.k2d',
    'gottcha-v-speDB': dirname(dirname(__file__)) + '/database/GOTTCHA/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.species',
    'gottcha-b-speDB': dirname(dirname(__file__)) + '/database/GOTTCHA/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species',
    'gottcha-v-strDB': dirname(dirname(__file__)) + '/database/GOTTCHA/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.strain',
    'gottcha-b-strDB': dirname(dirname(__file__)) + '/database/GOTTCHA/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.strain',
    'gottcha-v-genDB': dirname(dirname(__file__)) + '/database/GOTTCHA/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.genus',
    'gottcha-b-genDB': dirname(dirname(__file__)) + '/database/GOTTCHA/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.genus',
    'gottcha2-v-genDB': dirname(dirname(__file__)) + '/database/GOTTCHA2/RefSeq-Release89.Virus.genus.fna.gz',
    'gottcha2-b-speDB': dirname(dirname(__file__)) + '/database/GOTTCHA2/RefSeq-r90.cg.BacteriaArchaeaViruses.species.fna',
    'gottcha2-v-speDB': dirname(dirname(__file__)) + '/database/GOTTCHA2/RefSeq-Release90.cg.Viruses.species.fna',
    'gottcha2-e-plnDB': dirname(dirname(__file__)) + '/database/GOTTCHA2/RefSeq-Release89.Plant.species.fna.gz',
    'gottcha2-e-fugDB': dirname(dirname(__file__)) + '/database/GOTTCHA2/RefSeq-Release89.Fungi.species.fna.gz',
    'gottcha2-e-ptzDB': dirname(dirname(__file__)) + '/database/GOTTCHA2/RefSeq-Release89.Protozoa.species.fna.gz',
    'diamond-db': dirname(dirname(__file__)) + '/database/diamond/RefSeq_Release83.nr_protein_withRefSeq_viral_102317.protein.faa.dmnd',
    'centrifuge-db': dirname(dirname(__file__)) + '/database/Centrifuge/hpv.1.cf'
}

# PanGIA configs
config_json = json.loads(open(opt["configJson"]).read())
opt["pangia-db"] = config_json.get("edge-taxa-pangia-db", f"{EDGE_HOME}/database/PanGIA/NCBI_genomes_refseq89*.fa")
opt["pangia-bg"] = f"-b {EDGE_HOME}/database/PanGIA/background/{config_json['edge-taxa-pangia-bg']}" if config_json.get('edge-taxa-pangia-bg') else ""
opt["pangia-ra"] = config_json.get("edge-taxa-pangia-ra", "DEPTH_COV")
opt["pangia-ms"] = config_json.get("edge-taxa-pangia-ms", "0")
opt["pangia-mr"] = config_json.get("edge-taxa-pangia-mr", "3")
opt["pangia-mb"] = config_json.get("edge-taxa-pangia-mb", "1")
opt["pangia-ml"] = config_json.get("edge-taxa-pangia-ml", "50")
opt["pangia-rc"] = config_json.get("edge-taxa-pangia-rc", "R_MAT")
# opt["pangia-opts"] = ""
opt["pangia-opts"] = "-ps" if config_json.get("edge-taxa-pangia-ps-sw") else ""

tools = opt['tools'].split(',')
if not os.path.exists(opt['template']):
    print(usage(), file=sys.stderr)
if len(tools) < 1:
    print(usage(), file=sys.stderr)

# Remove suffix if any to meet the tool db format
for db in opt.keys():
    if db in ['bwa-db', 'gottcha.*DB', 'pangia']:
        opt[db] = re.sub(r'\.?(amb|ann|bwt|fai|pac|sa|parsedGOTTCHA\.dmp)$', '', opt[db])
    elif db == 'metaphlan-db':
        opt[db] = re.sub(r'(\.rev)?.\d\.bt2$', '', opt[db])
    elif db == 'centrifuge-db':
        opt[db] = re.sub(r'\.\d\.cf$', '', opt[db])
    elif db == 'kraken-db':
        tmp_filename, opt['kraken-db'], tmp_suffix = fileparse(opt['kraken-db'])

if opt['nanopore']:
    opt['gottcha-opts'] = "-a '-x ont2d'"
    opt['gottcha2-opts'] = "-a '--nanopore' "
    opt['bwa-opts'] = "-a '-x ont2d'"
    opt['pangia-opts'] = "-sb -se --nanopore "

for tool in tools:
    if tool:
        opt[tool] = 1

# print Dumper (\%opt)
template = HTMLTemplate(opt['template'], die_on_bad_params=0)
# template = HTMLTemplate(opt['template'], die_on_bad_params=1)
template.param(opt)
print(template.output())

def read_list_from_json(json_file):
    list_data = {}
    try:
        with open(json_file, 'r') as f:
            list_data = json.load(f)
    except FileNotFoundError:
        pass
    return list_data

def usage():
    print('''{} [template.tmpl] [tools] > microbial_profiling_configure.settings.ini 

    [Options]           [Default]
    -------------------------------------------------------------------
    -template           $Bin/microbial_profiling.settings.tmpl
    -tools              comma separated tools: bwa,kraken-mini, ...
    -bwaScoreCut        minimum score to output for BWA [30]
    -bwa-db             $EDGE_HOME/database/bwa_index/NCBI-Bacteria-Virus.fna
    -metaphlan-db       $EDGE_HOME/database/metaphlan/mpa
    -kraken-db          $EDGE_HOME/database/kraken2
    -centrifuge-db      $EDGE_HOME/database//Centrifuge/hpv.1.cf
    -gottcha-v-speDB    $EDGE_HOME/database/GOTTCHA/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.species
    -gottcha-b-speDB    $EDGE_HOME/database/GOTTCHA/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.species
    -gottcha-v-strDB    $EDGE_HOME/database/GOTTCHA/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.strain
    -gottcha-b-strDB    $EDGE_HOME/database/GOTTCHA/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.strain
    -gottcha-v-genDB    $EDGE_HOME/database/GOTTCHA/GOTTCHA_VIRUSES_c5900_k24_u30_xHUMAN3x.genus
    -gottcha-b-genDB    $EDGE_HOME/database/GOTTCHA/GOTTCHA_BACTERIA_c4937_k24_u30_xHUMAN3x.genus
    -gottcha2-v-genDB   $EDGE_HOME/database/GOTTCHA2/RefSeq-Release89.Virus.genus.fna.gz
    -gottcha2-b-speDB   $EDGE_HOME/database/GOTTCHA2/RefSeq-r90.cg.BacteriaViruses.species.fna
    -gottcha2-v-speDB   $EDGE_HOME/database/GOTTCHA2/RefSeq-Release89.Virus.species.fna.gz
    -gottcha2-e-plnDB   $EDGE_HOME/database/GOTTCHA2/RefSeq-Release89.Plant.species.fna.gz
    -gottcha2-e-ptzDB   $EDGE_HOME/database/GOTTCHA2/RefSeq-Release89.Protozoa.species.fna.gz
    -gottcha2-e-fugDB   $EDGE_HOME/database/GOTTCHA2/RefSeq-Release89.Fungi.species.fna.gz
    -diamond-db         $EDGE_HOME/database/diamond/RefSeq_Release83.nr_protein_withRefSeq_viral_102317.protein.faa.dmnd
    --nanopore          

'''.format(sys.argv[0]))
    sys.exit(0)

if __name__ == '__main__':
    usage()
    sys.exit(0)
