#!/usr/bin/env python3
# sra2fastq.py
# ver 0.5
# 2014/12/19
#
# Po-E (Paul) Li
# B-11
# Los Alamos National Lab.

# Change log
# ver 0.5 (2016/11/02)
# - support downloading sequences from NCBI-SRA, EBI-ENA and DDBJ.
# - auto switch download sites.
# ver 0.3 (2016/10/28)
# - Switch sequence source to ERA
# ver 0.2 (2015/01/06)
# - Input SRA accessions support studies (SRP*/ERP*/DRP*), experiments (SRX*/ERX*/DRX*), samples (SRS*/ERS*/DRS*), runs (SRR*/ERR*/DRR*), or submissions (SRA*/ERA*/DRA*)
# - Provide "--platform-restrict" option to limit the platform of SRAs
# - Provide "--concat" option to concatenate multiple FASTQ files into a singal (single-end) or two (paired-end) files
# - Remove dependency of File::Which
# ver 0.1
# - Initial release

import os
import re
import subprocess
import sys
import shutil
from getopt import getopt
from pathlib import Path


OUTDIR = '.'
platform_restrict = None
clean = False
filesize_restrict = 0
runs_restrict = 0
download_tool = 'curl'
user_proxy = None
no_proxy = False
http_proxy = os.environ.get('HTTP_PROXY') or os.environ.get('http_proxy')
ftp_proxy = os.environ.get('FTP_PROXY') or os.environ.get('ftp_proxy')

http_proxy = f"--proxy '{http_proxy}' " if http_proxy else ""
ftp_proxy = f"--proxy '{ftp_proxy}' " if ftp_proxy else ""
gzip = None



## Subroutines ########################################################################

def getSraFastq(info, run_acc):
    print(f"Retrieving FASTQ for {run_acc} from NCBI SRA (online converting)...")
    platform = info['platform']
    library = info['library']
    url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra-reads-be/fastq?acc={run_acc}"
    cmd = f"{curl} -o {OUTDIR}/sra2fastq_temp/{run_acc}.fastq.gz \"{url}\""
    
    print(f"Downloading {url}...")
    ec = subprocess.call(cmd, shell=True)
    
    if ec > 0:
        print(f"Failed to download SRA file from {url}.")
        return "failed"
    
    # Deinterleaving if paired-end reads
    if "illu" in platform.lower() and "pair" in library.lower():
        print("Paired-end reads found. Deinterleaving...")
        di_flag = subprocess.call(f"gzip -dc {OUTDIR}/sra2fastq_temp/{run_acc}.fastq.gz | deinterleave_fastq.sh {OUTDIR}/sra2fastq_temp/{run_acc}_1.fastq.gz {OUTDIR}/sra2fastq_temp/{run_acc}_2.fastq.gz compress", shell=True)
        if di_flag > 0:
            return "failed"
        print("Done.")
        subprocess.call(f"rm -f {OUTDIR}/sra2fastq_temp/{run_acc}.fastq.gz", shell=True)
    
    return "success"



def getSraFastqToolkits(info, run_acc):
    print(f"Retrieving FASTQ for {run_acc} with NCBI SRA Toolkit...")

    platform = info["platform"]
    url = info["url"]
    filename = run_acc

    print(f"Downloading {url}...")
    cmd = f"{curl} -O {OUTDIR}/sra2fastq_temp/{filename} \"{url}\"" if download_tool == "wget" else f"{curl} {http_proxy} -o {OUTDIR}/sra2fastq_temp/{filename} \"{url}\""
    ec = subprocess.call(cmd, shell=True)
    
    if ec > 0:
        print(f"Failed to download SRA file from {url}.")
        return "failed"
    print("Done.")

    #check downloaded file
    filesize = os.path.getsize(f"{OUTDIR}/sra2fastq_temp/{filename}")

    if not filesize:
        print(f"Failed to download SRA file from {url}.")
        return "failed"

    #dump fastq from SRA file
    options = "--gzip "
    options += "--split-files " if "illu" in platform.lower() else ""
    options += "--split-files -B " if "solid" in platform.lower() else ""
    print(f"Running fastq-dump with options {options}...")
    ec = subprocess.call(f"fastq-dump {options} --outdir '{OUTDIR}/sra2fastq_temp' {OUTDIR}/sra2fastq_temp/{filename} 2>/dev/null", shell=True)
    
    if ec > 0:
        print(f"Failed to run fastq-dump from {OUTDIR}/sra2fastq_temp/{filename}.")
        return "failed"
    print("Done.")

    #clean up temp files
    os.remove(f"{OUTDIR}/sra2fastq_temp/{filename}")
    return "success"


def getDdbjFastq(info, run_acc):
    print("Retrieving FASTQ for", run_acc, "from DDBJ...", file=sys.stderr)

    platform = info['platform']
    library = info['library']
    exp_acc = info['exp_acc']
    sub_acc = info['sub_acc']
    sra_acc_first6 = sub_acc[:6]
    ec = None
    cmd = None
    MinSize = 10000

    url = f"ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/{sra_acc_first6}/{sub_acc}/{exp_acc}"

    if "illu" in platform.lower() and "pair" in library.lower():
        print(f"Downloading {url}/{run_acc}_1.fastq.bz2...", file=sys.stderr)
        cmd = f"{download_tool} -O {OUTDIR}/sra2fastq_temp/{run_acc}_1.fastq.bz2 \"{url}/{run_acc}_1.fastq.bz2\"" if download_tool == "wget" else f"{curl} {ftp_proxy} -o {OUTDIR}/sra2fastq_temp/{run_acc}_1.fastq.bz2 \"{url}/{run_acc}_1.fastq.bz2\""
        ec = subprocess.call(cmd, shell=True)
        print("finished.", file=sys.stderr)
        print(f"Downloading {url}/{run_acc}_2.fastq.bz2...", file=sys.stderr)
        cmd = f"{download_tool} -O {OUTDIR}/sra2fastq_temp/{run_acc}_2.fastq.bz2 \"{url}/{run_acc}_2.fastq.bz2\"" if download_tool == "wget" else f"{curl} {ftp_proxy} -o {OUTDIR}/sra2fastq_temp/{run_acc}_2.fastq.bz2 \"{url}/{run_acc}_2.fastq.bz2\""
        ec = subprocess.call(cmd, shell=True)
        print("finished.", file=sys.stderr)

    print(f"Downloading {url}/{run_acc}.fastq.bz2...", file=sys.stderr)
    cmd = f"{curl} {ftp_proxy} -o {OUTDIR}/sra2fastq_temp/{run_acc}.fastq.bz2 \"{url}/{run_acc}.fastq.bz2\""
    ec = subprocess.call(cmd, shell=True)
    print("finished.", file=sys.stderr)

    total_size = 0
    if os.path.getsize(f"{OUTDIR}/sra2fastq_temp/{run_acc}_1.fastq.bz2") > MinSize:
        print("Convering bz2 to gz...")
        ec = os.system(f"bunzip2 -c < {OUTDIR}/sra2fastq_temp/{run_acc}_1.fastq.bz2 | gzip -c > {OUTDIR}/sra2fastq_temp/{run_acc}_1.fastq.gz")
        if ec > 0:
            print("failed to convert bz2 to gz.")
            return "failed"
        total_size += os.path.getsize(f"{OUTDIR}/sra2fastq_temp/{run_acc}_1.fastq.gz")
        print("Done.")
    else:
        os.unlink(f"{OUTDIR}/sra2fastq_temp/{run_acc}_1.fastq.bz2")

    if os.path.getsize(f"{OUTDIR}/sra2fastq_temp/{run_acc}_2.fastq.bz2") > MinSize:
        print("Convering bz2 to gz...")
        ec = os.system(f"bunzip2 -c < {OUTDIR}/sra2fastq_temp/{run_acc}_2.fastq.bz2 | gzip -c > {OUTDIR}/sra2fastq_temp/{run_acc}_2.fastq.gz")
        if ec > 0:
            print("failed to convert bz2 to gz.")
            return "failed"
        total_size += os.path.getsize(f"{OUTDIR}/sra2fastq_temp/{run_acc}_2.fastq.gz")
        print("Done.")
    else:
        os.unlink(f"{OUTDIR}/sra2fastq_temp/{run_acc}_2.fastq.bz2")

    if os.path.getsize(f"{OUTDIR}/sra2fastq_temp/{run_acc}.fastq.bz2") > MinSize:
        print("Convering bz2 to gz...")
        ec = os.system(f"bunzip2 -c < {OUTDIR}/sra2fastq_temp/{run_acc}.fastq.bz2 | gzip -c > {OUTDIR}/sra2fastq_temp/{run_acc}.fastq.gz")
        if ec > 0:
            print("failed to convert bz2 to gz.")
            return "failed"
        total_size += os.path.getsize(f"{OUTDIR}/sra2fastq_temp/{run_acc}.fastq.gz")
        print("Done.")
    else:
        os.unlink(f"{OUTDIR}/sra2fastq_temp/{run_acc}.fastq.bz2")

    print(total_size)

    if total_size < 50:
        os.system(f"rm {OUTDIR}/sra2fastq_temp/{run_acc}*gz")
        print("failed to download FASTQ and convert FASTQ files from DDBJ.")
        return "failed"

    return "success"


def getEnaFastq(info, run_acc):
    print("Retrieving FASTQ for", run_acc, "from EBI-ENA...")
    for i in info['ena_fastq_ftp']:
        url = info['ena_fastq_ftp'][i]['url']
        md5 = info['ena_fastq_ftp'][i]['md5']
        size = info['ena_fastq_ftp'][i]['size']
        
        filename = url.split('/')[-1]
        print("Downloading", url, "...")
        
        # using subprocess.call to execute curl command
        cmd = f"{curl} -o {os.path.join(OUTDIR, 'sra2fastq_temp', filename)} \"ftp://{url}\""
        subprocess.call(cmd, shell=True)
        
        # check downloaded file
        filesize = os.path.getsize(os.path.join(OUTDIR, 'sra2fastq_temp', filename))
        if not filesize:
            print(f"Failed to download {filename} from {url}.")
            return "failed"
        if filesize != size:
            print(f"{os.path.join(OUTDIR, 'sra2fastq_temp', filename)} incomplete/corrupted -- file sizes mismatch.")
            return "failed"
        
        # check md5
        output = subprocess.check_output(f"md5sum {os.path.join(OUTDIR, 'sra2fastq_temp', filename)}", shell=True)
        md5sum = output.split()[0].decode('utf-8')
        if md5sum != md5:
            print(f"{os.path.join(OUTDIR, 'sra2fastq_temp', filename)} corrupted -- md5 checksum mismatch.")
            return "failed"
        
        print("Done.")
    
    return "success"



def get_read_info(acc, read_info, sra_type):
    print("Retrieving run(s) information from NCBI-SRA...")
    
    # get info from NCBI-SRA
    url0 = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={acc}&usehistory=y"
    download_tool = 'wget' # or any other download tool
    cmd = f"{download_tool} -O - \"{url0}\" 2>/dev/null" if download_tool == 'wget' else f"{download_tool} {http_proxy} \"{url0}\" 2>/dev/null"
    web_result = subprocess.check_output(cmd, shell=True).decode('utf-8')
    lines = web_result.split("\n")
    webenv = ''
    key = ''
    for line in lines:
        if re.search('<WebEnv>(\S+)<\/WebEnv>', line):
            webenv = re.search('<WebEnv>(\S+)<\/WebEnv>', line).group(1)
        if re.search('<QueryKey>(\S+)<\/QueryKey>', line):
            key = re.search('<QueryKey>(\S+)<\/QueryKey>', line).group(1)
    
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&rettype=runinfo&query_key={key}&WebEnv={webenv}&retmode=text"
    print(f"Retrieving run acc# from NCBI-SRA {url0} {url}...")
    cmd = f"{download_tool} -O - \"{url}\" 2>/dev/null" if download_tool == 'wget' else f"{download_tool} {http_proxy} \"{url}\" 2>/dev/null"
    web_result = subprocess.check_output(cmd, shell=True).decode('utf-8')
    lines = web_result.split("\n")
    sra_num_runs = len(lines)
    print(f"{sra_num_runs} run(s) found from NCBI-SRA.")

    for line in lines:
        if line.startswith("Run"):
            continue
        fields = line.split(',')
        sub_acc = fields[42] #submission
        exp_acc = fields[10] #experiment_accession
        run_acc = fields[0] #run
        size_MB = fields[7]
        platform = fields[18] #platform
        library = fields[15] #LibraryLayout
        url = fields[9] #download_path
        if sra_type == "ByRun" and not re.search(acc, run_acc, re.IGNORECASE):
            continue
        if not size_MB:
            print(f"Run {run_acc} has size 0 MB")
        read_info[acc][run_acc] = {
            "exp_acc": exp_acc,
            "sub_acc": sub_acc,
            "platform": platform,
            "library": library,
            "url": url
        }

    ##### get info from EBI-ENA when NCBI-SRA fails ####

    print("Retrieving run(s) information from EBI-ENA...\n", file=sys.stderr)

    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run&fields=run_accession,submission_accession,study_accession,experiment_accession,instrument_platform,library_layout,fastq_ftp,fastq_md5,fastq_bytes"
    print(f"Retrieving run acc# from EBI-ENA {url}...\n", file=sys.stderr)
    cmd = f"{download_tool} {http_proxy} {url}" if download_tool != "wget" else f"wget -O - {url}"
    web_result = subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL).decode()

    if not re.search(r"^study_accession|^run_accession", web_result) and not sra_num_runs:
        raise Exception(f"ERROR: Failed to retrieve sequence information for {acc} from both SRA and ENA database.")
    elif not re.search(r"^study_accession|^run_accession", web_result) and sra_num_runs:
        print(f"WARNING: {acc} only found in SRA database. The data may be not synchronized among INSDC yet.\n", file=sys.stderr)
    elif re.search(r"^study_accession|^run_accession", web_result) and not sra_num_runs:
        print(f"WARNING: {acc} only found in ENA database. The data may be not synchronized among INSDC yet.\n", file=sys.stderr)

    lines = web_result.split('\n')
    print(f"{len(lines)-2} run(s) found from EBI-ENA.\n", file=sys.stderr)

    for line in lines:
        if re.search(r"^study_accession|^run_accession", line):
            continue
        fields = line.strip().split('\t')

        sub_acc = fields[1] # submission_accession
        exp_acc = fields[3] # experiment_accession
        run_acc = fields[0] # run_accession
        platform = fields[4] # instrument_platform
        library = fields[5] # library_layout

        url = fields[6].split(';') # fastq_ftp
        md5 = fields[7].split(';') # fastq_md5
        size = fields[8].split(';') # fastq_bytes

        if sra_type == "ByRun" and not re.search(acc, run_acc, re.IGNORECASE):
            continue

        read_info.setdefault(acc, {}).setdefault(run_acc, {})["platform"] = platform
        read_info.setdefault(acc, {}).setdefault(run_acc, {})["exp_acc"] = exp_acc
        read_info.setdefault(acc, {}).setdefault(run_acc, {})["sub_acc"] = sub_acc
        read_info.setdefault(acc, {}).setdefault(run_acc, {})["library"] = library

        for i in range(len(url)):
            read_info.setdefault(acc, {}).setdefault(run_acc, {}).setdefault("fastq_ftp", {}).setdefault(i, {})["url"] = url[i]
            read_info.setdefault(acc, {}).setdefault(run_acc, {}).setdefault("fastq_ftp", {}).setdefault(i, {})["md5"] = md5[i]
            read_info.setdefault(acc, {}).setdefault(run_acc, {}).setdefault("fastq_ftp", {}).setdefault(i, {})["size"] = size[i]
    
    return read_info

def check_required_exec():
    if subprocess.call(['which', 'gzip']) != 0:
        sys.exit("ERROR: 'gzip' not found.")
    if download_tool == 'curl' and subprocess.call(['which', 'curl']) != 0:
        sys.exit("ERROR: 'curl' not found.")
    if download_tool == 'wget' and subprocess.call(['which', 'wget']) != 0:
        sys.exit("ERROR: 'wget' not found.")

def usage():
    print('''
[DESCRIPTION]
    A script retrieves sequence project in FASTQ files from 
NCBI-SRA/EBI-ENA/DDBJ database using `curl` or `wget`. Input accession number
supports studies (SRP*/ERP*/DRP*), experiments (SRX*/ERX*/DRX*), 
samples (SRS*/ERS*/DRS*), runs (SRR*/ERR*/DRR*), or submissions 
(SRA*/ERA*/DRA*).

[USAGE]
    $0 [OPTIONS] <Accession#> (<Accession# 2> <Accession# 3>...)

[OPTIONS]
    --outdir|-d            Output directory
    --clean                Clean up temp directory
    --platform-restrict    Only allow a specific platform
    --filesize-restrict    (in MB) Only allow to download less than a specific
                           total size of files.
    --run-restrict         Only allow download less than a specific number
                           of runs.
    --download-interface   curl or wget [default: curl]
    --help/-h/-?           Display this help
''')
    sys.exit()

## Check Env ########################################################################

check_required_exec()

opts, args = getopt(
    sys.argv[1:],
    'd:pr:fr:r:n',
    [
        'OUTDIR=',
        'platform-restrict=',
        'filesize-restrict=',
        'runs-restrict=',
        'download-interface=',
        'proxy=',
        'no_proxy',
        'clean',
        'help'
    ]
)

for opt, val in opts:
    if opt in ('-d', '--OUTDIR'):
        OUTDIR = val
    elif opt in ('-pr', '--platform-restrict'):
        platform_restrict = val
    elif opt in ('-fr', '--filesize-restrict'):
        filesize_restrict = val
    elif opt in ('-r', '--runs-restrict', '-rr'):
        runs_restrict = val
    elif opt == '--download-interface':
        download_tool = val
    elif opt == '--proxy':
        user_proxy = val
    elif opt == '--no_proxy':
        no_proxy = True
    elif opt in ('-n', '--clean'):
        clean = True
    elif opt == '--help':
        usage()
        sys.exit()

if not args:
    usage()
    sys.exit()

http_proxy = "" if no_proxy else http_proxy
ftp_proxy = "" if no_proxy else ftp_proxy
http_proxy = f"--proxy '{user_proxy}' " if user_proxy else http_proxy
ftp_proxy = f"--proxy '{user_proxy}' " if user_proxy else ftp_proxy

curl = f"wget -v -U 'Mozilla/5.0' " if download_tool == 'wget' else f"curl -k -A 'Mozilla/5.0' -L "

# init temp directory
subprocess.run(f'rm -rf {OUTDIR}/sra2fastq_temp/', shell=True, check=True)
Path(f'{OUTDIR}/sra2fastq_temp/merged/').mkdir(parents=True, exist_ok=True)




## Main ########################################################################

read_info = {}
total_runs = 0
total_size = 0
dl_status = None
dl_runs = None
dl_size = None

for acc in sys.argv[1:]:
    if not re.match(r'^(SR|ER|DR)', acc):
        sys.exit(f"ERROR: {acc} is not a valid SRA/ERA/DRA number.")

    sra_type = "ByRun"
    if re.match(r'^(SRX|ERX|DRX)', acc):
        sra_type = "ByExp"
    elif re.match(r'^(SRS|ERS|DRS)', acc):
        sra_type = "BySample"
    elif re.match(r'^(SRP|ERP|DRP)', acc):
        sra_type = "ByStudy"

    print(f"Processing {acc} ({sra_type})...")

    # gathering reads information from NCBI-SRA / EBI-ENA
    read_info = get_read_info(acc, read_info, sra_type)

    if acc not in read_info:
        sys.exit(f"ERROR: No sequence found. Please check if {acc} is a valid SRA/ERA/DRA number or your internet connection.")

    for run_acc in read_info[acc]:
        # check if runs exceed the limit
        total_runs += 1
        if runs_restrict and runs_restrict < total_runs:
            sys.exit("ERROR: Run(s) exceed the limit ({} MB).".format(runs_restrict))

        # check platform
        platform = read_info[acc][run_acc]['platform']
        if platform_restrict and not re.search(platform_restrict, platform, re.IGNORECASE):
            print(f"WARN: {read_info[acc][run_acc]['platform']} platform detected. Only {platform_restrict} is allowed.")
            continue

        # download FASTQ
        dl_status = getSraFastqToolkits(read_info[acc][run_acc], run_acc)
        if dl_status == 'failed':
            dl_status = getSraFastq(read_info[acc][run_acc], run_acc)
        if dl_status == 'failed':
            dl_status = getDdbjFastq(read_info[acc][run_acc], run_acc)
        if dl_status == 'failed':
            dl_status = getEnaFastq(read_info[acc][run_acc], run_acc)
        if dl_status == 'failed':
            sys.exit("ERROR: Please check your internet connection.")



        # merging fastqs in multiple runs into a single file
        if os.path.exists(f"{OUTDIR}/sra2fastq_temp/{run_acc}_1.fastq.gz"):
            with open(f"{OUTDIR}/sra2fastq_temp/merged/{acc}.1.fastq.gz", "ab") as out_file:
                with open(f"{OUTDIR}/sra2fastq_temp/{run_acc}_1.fastq.gz", "rb") as in_file:
                    out_file.write(in_file.read())
                total_size += os.path.getsize(f"{OUTDIR}/sra2fastq_temp/{run_acc}_1.fastq.gz")
            os.remove(f"{run_acc}_1.fastq.gz")

        if os.path.exists(f"{OUTDIR}/sra2fastq_temp/{run_acc}_2.fastq.gz"):
            with open(f"{OUTDIR}/sra2fastq_temp/merged/{acc}.2.fastq.gz", "ab") as out_file:
                with open(f"{OUTDIR}/sra2fastq_temp/{run_acc}_2.fastq.gz", "rb") as in_file:
                    out_file.write(in_file.read())
                total_size += os.path.getsize(f"{OUTDIR}/sra2fastq_temp/{run_acc}_2.fastq.gz")
            os.remove(f"{run_acc}_2.fastq.gz")

        if os.path.exists(f"{OUTDIR}/sra2fastq_temp/{run_acc}.fastq.gz"):
            with open(f"{OUTDIR}/sra2fastq_temp/merged/{acc}.fastq.gz", "ab") as out_file:
                with open(f"{OUTDIR}/sra2fastq_temp/{run_acc}.fastq.gz", "rb") as in_file:
                    out_file.write(in_file.read())
                total_size += os.path.getsize(f"{OUTDIR}/sra2fastq_temp/{run_acc}.fastq.gz")
            os.remove(f"{run_acc}.fastq.gz")

        if filesize_restrict and filesize_restrict < total_size/1024/1024:
            raise ValueError("ERROR: downloaded file size exceed limitation ({filesize_restrict} MB).\n")
        else:
            print("Succesfully downloaded {run_acc}.")
    
    sys.stderr.write(f"Finished downloading acc# {acc}.\n")

subprocess.run(["mv", f"{OUTDIR}/sra2fastq_temp/merged/*.gz", OUTDIR], check=True)

if clean:
    print("cleaning up...")
    shutil.rmtree(os.path.join(OUTDIR, "sra2fastq_temp"))
    print("Done.")