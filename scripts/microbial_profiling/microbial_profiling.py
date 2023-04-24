#!/usr/bin/env python
# Paul Li, B-11 (April 2013)
#

import os
import sys
import getopt
import subprocess
import shlex

# Constants
threads = 4
extract = 0.1

# Define usage function
def usage():
    print("Usage: ")
    print("\tperl2python.pl [options]")
    print("Options:")
    print("\t--list, -l: input file list")
    print("\t--setting, -s: settings file")
    print("\t--cpu, -c: number of CPU threads (default: 4)")
    print("\t--output, -o: output directory (default: current directory)")
    print("\t--debug, -d: enable debug mode")
    print("\t--help, -h, -?: display help message")

# Parse command line options
try:
    opts, args = getopt.getopt(sys.argv[1:], "l:s:c:o:dh?", ["list=", "setting=", "cpu=", "output=", "debug", "help"])
except getopt.GetoptError as err:
    print(str(err))
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ("-l", "--list"):
        list_file = arg
    elif opt in ("-s", "--setting"):
        setting_file = arg
    elif opt in ("-c", "--cpu"):
        threads = int(arg)
    elif opt in ("-o", "--output"):
        output_dir = arg
    elif opt in ("-d", "--debug"):
        debug = True
    elif opt in ("-h", "--help", "-?"):
        usage()
        sys.exit()

# Check for required options
if not 'setting_file' in locals():
    print("ERROR: No settings file provided.")
    usage()
    sys.exit(2)

# Read settings
ini_file = setting_file
tools = restore_settings(ini_file)
pid = os.getpid()

# Default settings
threads = tools['system']['THREADS'] if 'THREADS' in tools['system'] else threads
if 'cpu' in locals():
    threads = cpu
extract = tools['system']['EXTRACT_NUM'] if 'EXTRACT_NUM' in tools['system'] else extract

#sequedex out script
sqdx_fh = None
post_fh = open(post_script, "w") or sys.exit("ERROR: Can't create post-processing script file {}: {}".format(post_script, os.strerror(errno.EEXIST)))
info_fh = None
if not os.path.exists(fileinfo_out):
    info_fh = open(fileinfo_out, "w") or sys.exit("ERROR: Can't create fileinfo file {}: {}".format(fileinfo_out, os.strerror(errno.EEXIST)))
list_fh = None
if not os.path.exists(opt['list']):
    list_fh = open(filelist_out, "w") or sys.exit("ERROR: Can't create filelist {}: {}".format(filelist_out, os.strerror(errno.EEXIST)))

log_fh = None
if logfile is not None:
    opt['verbose'] = 1
    log_fh = open(logfile, "w") or sys.exit("ERROR: Can't create log file {}: {}".format(logfile, os.strerror(errno.EEXIST)))

# retrieve input files
file_info = {}
cmd = ""
count = 1

files = argv

filelist_header = []
if opt['list'] is not None:
    with open(opt['list'], 'r') as LIST:
        for line in LIST:
            line = line.strip()
            if line.startswith("--") or line.strip() == "":
                continue
            if line.startswith("PREFIX"):
                filelist_header = line.split("\t")
                continue
            if len(filelist_header) > 0:
                fields = line.split("\t")
                for i in range(len(filelist_header)):
                    file_info[count][filelist_header[i]] = fields[i]
            else:
                files.append(line)
            count += 1

for file in files:
    if "," in file:
        file_info[count]["FASTQPE"] = file
    else:
        file_info[count]["FASTQ"] = file
    count += 1

num = len(file_info)
if num < 1:
    raise ValueError("No input FASTQ files.")

prepSequence(file_info, tools)

outbase = os.path.basename(files[0])
outbase = outbase.replace(".fq", "").replace(".fastq", "").replace(".gz", "").replace(".bz2", "")

post_script_content = "#!/bin/bash\n"
post_script_content += "module purge\n"
post_script_content += "module load perl\n"
post_script_content += "module load R/3.4.0\n"
post_script_content += "module load python/2.7.13\n"
post_script_content += "export LD_LIBRARY_PATH=/cm/shared/apps/gcc/4.8.2/lib64:$LD_LIBRARY_PATH\n"
post_script_content += "cd {}\n".format(p_outdir)
post_script_content += "R CMD BATCH --vanilla --slave --no-save --no-restore \"--args"
post_script_content += " --in={}".format(fileinfo_out)
post_script_content += " --prefix={}".format(outbase)
post_script_content += " --log={}".format(summary_out)
post_script_content += " --tool_list={}".format(tools['system']['SQDX_TOOL'])
post_script_content += " --report_template={}".format(tools['system']['REPORT_TEMPLATE'])
post_script_content += " --output_dir={}".format(p_repdir)
post_script_content += " --max_process_num={}".format(max_process_num)
post_script_content += "  # Add more arguments here if needed\"\n"

post_fh.write(post_script_content)
post_fh.close()

# run tools
_notify("\n[TOOLS]\n\n")
childs = []
job_ids = {}
toolnames = sorted(tools.keys(), key=lambda x: tools[x]['ORDER'])
if tools['system']['RUN_TOOLS']:
    for idx in sorted(file_info.keys()):
        for tool in toolnames:
            if tool == 'system':
                continue
            input = file_info[idx]['FASTQ']
            fnb = file_info[idx]['PREFIX']
            outdir = f"{p_outdir}/{idx}_{fnb}/{tool}"
            tool_rep_dir = f"{p_repdir}/{idx}_{fnb}/{tool}"
            prefix = fnb
            log = f"{p_logdir}/{fnb}-{tool}.log"
            print(f"Tool ({tool}) - PID: {os.getpid()}, starting...")

            # skip if tool output exists
            if os.path.exists(f"{tool_rep_dir}/.finished"):
                _notify(f"[RUN_TOOL] [{tool}] Result exists. Skipping tool {tool}!\n")
                _notify(f"[RUN_TOOL] [{tool}] Running time: 00:00:00\n")
                continue

            # prepare command
            code = None
            time = int(time.time())
            cmd = tools[tool]['COMMAND']
            qsub_cmd = tools['system']['QSUB_COMMAND']
            cmd = param_replace(cmd, file_info, tools, idx, tool)
            qsub_cmd = param_replace(qsub_cmd, file_info, tools, idx, tool)
            if tools['system']['RUN_TOOL_AS_JOB']:
                os.unlink(log)
                _notify(f"[RUN_TOOL] [{tool}] COMMAND: {qsub_cmd}\n")
                _notify(f"[RUN_TOOL] [{tool}] Logfile: {log}\n")
                job_id = None
                if 'qsub' in qsub_cmd:
                    code = subprocess.check_output(f"{qsub_cmd} -v EDGE_HOME=$ENV['EDGE_HOME']} {script_dirname}/script/{cmd}", shell=True, stderr=subprocess.STDOUT, universal_newlines=True)
                    match = re.search(r'Your job (\d+)', code)
                    if match:
                        job_id = match.group(1)
                elif 'sbatch' in qsub_cmd:
                    code = subprocess.check_output(f"{qsub_cmd} --export=EDGE_HOME=$ENV['EDGE_HOME']} --wrap='{script_dirname}/script/{cmd}'", shell=True, stderr=subprocess.STDOUT, universal_newlines=True)
                    match = re.search(r'Submitted batch job (\d+)', code)
                    if match:
                        job_id = match.group(1)
                _notify(f"[RUN_TOOL] [{tool}] cluster info: {code}\n")
                if job_id:
                    job_ids[job_id] = {'tool': tool, 'time': time, 'outdir': outdir, 'log': log}
                else:
                    _notify(f"[RUN_TOOL] [{tool}] Error occurred: qsub error\n")
            else:
                _notify(f"[RUN_TOOL] [{tool}] COMMAND: {cmd}\n")
                _notify(f"[RUN_TOOL] [{tool}] Logfile: {log}\n")
                code = os.system(f"{cmd} > {log} 2>&1")
                if code:
                    if 'No target reads' in open(log).read():
                        _notify(f"[RUN_TOOL] [{tool}] Warning: No target reads found.\n")
                    else:
                        _notify(f"[RUN_TOOL] [{tool}] Error occurred.\n")
                runningtime = timeInterval(time)
                _notify(f"[RUN_TOOL] [{tool}] Running time: {runningtime}\n")
                if not code:
                    os.system(f"touch \"{outdir}/.finished\"")
else:
    _notify("[RUN_TOOL] Skipped.\n")


if tools['system']['RUN_TOOL_AS_JOB']:
    timelimit = tools['system']['QSUB_TIMELIMIT']
    wait = True
    qstat_cmd = tools['system']['QSTAT_COMMAND']
    qdel_cmd = tools['system']['QDEL_COMMAND']
    while wait:
        qsub_job_exit = 0
        for id in job_ids.keys():
            if job_ids[id]['finish']:
                qsub_job_exit += 1
                continue
            job_status = subprocess.check_output([qstat_cmd, id], stderr=subprocess.STDOUT).decode('utf-8')
            time = job_ids[id]['time']
            tool = job_ids[id]['tool']
            outdir = job_ids[id]['outdir']
            log = job_ids[id]['log']
            if 'do not exist' in job_status or 'Invalid job id specified' in job_status:
                qsub_job_exit += 1
                runningtime = timeInterval(time)
                if 'No target reads' in open(log).read():
                    notify("[RUN_TOOL] [{}] Warning: No target reads found.\n".format(tool))
                elif 'ERROR' in open(log).read():
                    notify("[RUN_TOOL] [{}] Error occurred.\n".format(tool))
                else:
                    open(outdir + '/.finished', 'a').close()
                notify("[RUN_TOOL] [{}] Running time: {}\n".format(tool, runningtime))
                job_ids[id]['finish'] = 1
            if (time.time() - time) > timelimit:
                subprocess.call([qdel_cmd, id])
                qsub_job_exit += 1
                job_ids[id]['finish'] = 1
                notify("[RUN_TOOL] [{}] Error: TIME OUT\n".format(tool))
        if len(job_ids.keys()) == qsub_job_exit:
            break
        if not job_ids:
            break
        time.sleep(5)

for child in childs:
    os.waitpid(child, 0)


# generate post-process script
heatmap_scale = tools.get('system', {}).get('HEATMAP_SCALE', 'log')
heatmap_top = tools.get('system', {}).get('HEATMAP_DISPLAY_TOP', 0)

for idx in sorted(file_info.keys()):
    fa = file_info[idx]['FASTA']
    fnb = file_info[idx]['PREFIX']
    tmpdir = f"{p_outdir}/temp"
    gottcha_present = 0

    pwd = os.getcwd()

    print(f"""
      export PATH={Bin}:{Bin}/../:{Bin}/script/:{Bin}/bin/:{os.environ['EDGE_HOME']}/thirdParty/Anaconda2/bin:$PATH;
      cd {pwd};

      echo "[Post-processing #{idx} {fnb}]";
      mkdir -p {tmpdir};
	""")

    for tool in sorted(tools.keys(), key=lambda x: tools[x]['ORDER']):
        if tool == 'system':
            continue

        print("\n(\n")

        if 'gottcha' in tool.lower():
            gottcha_present = 1

        outdir = f"{p_outdir}/{idx}_{fnb}/{tool}"
        tool_rep_dir = f"{p_repdir}/{idx}_{fnb}/{tool}"
        prefix = f"{fnb}"

        print(f"echo \"==> processing result: {tool}\";\n")

        # copy output list & krona files
        post_fh.write(f"""
          mkdir -p {tool_rep_dir}
          echo "====> Copying result list to report directory..."
          if [ -e "{outdir}/{prefix}.out.list" ]; then
              cp {outdir}/{prefix}.out.list {tool_rep_dir}/{fnb}-{tool}.list.txt
          fi
          if [ -e "{outdir}/{prefix}.out.megan" ]; then
              cp {outdir}/{prefix}.out.megan {tool_rep_dir}/{fnb}-{tool}.megan
          fi
          if [ -e "{outdir}/{prefix}.full.tsv" ]; then
              cp {outdir}/{prefix}.full.tsv {tool_rep_dir}/{fnb}-{tool}.full.tsv
          fi
          if [ -e "{outdir}/{prefix}.krona.html" ]; then
              cp {outdir}/{prefix}.krona.html {tool_rep_dir}/{fnb}-{tool}.krona.html
          fi
          if [ -e "{outdir}/{prefix}.classification.csv" ]; then
              cp {outdir}/{prefix}.classification.csv {tool_rep_dir}/{fnb}-{tool}.classification.csv
          fi
          if [ -e "{outdir}/.finished" ]; then
              cp {outdir}/.finished {tool_rep_dir}/.finished
          fi
          if [ -e "{outdir}/{prefix}.sam" ]; then
              subprocess.run(["samtools", "view", "-b", "-@", threads, "-S", f"{outdir}/{prefix}.sam", "-o", f"{tool_rep_dir}/{fnb}-{tool}.bam"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
          fi
          if [ -e "{outdir}/{prefix}.gottcha.sam" ]; then
              subprocess.run(["samtools", "view", "-b", "-@", threads, "-S", f"{outdir}/{prefix}.gottcha.sam", "-o", f"{tool_rep_dir}/{fnb}-{tool}.bam"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
          fi
          if [ -e "{outdir}/{prefix}.pangia.sam" ]; then
              shutil.copytree(outdir, tool_rep_dir)
          fi
          if ls {outdir}/{prefix}.gottcha_*.sam 1>/dev/null 2>&1; then
              cp -f {outdir}/{prefix}.gottcha_*.sam {tool_rep_dir}/
          fi
          if [ -e "{outdir}/{prefix}.out.read_classification" ]; then
              cp {outdir}/{prefix}.out.read_classification {tool_rep_dir}/{fnb}-{tool}.read_classification
          fi
          echo "====> Generating phylo_dot_plot for each tool..."
          if [ -e "{outdir}/{prefix}.out.tab_tree" ]; then
              subprocess.run(["phylo_dot_plot.pl", "-i", f"{outdir}/{prefix}.out.tab_tree", "--score", f"{outdir}/{prefix}.out.tab_tree.score", "-p", f"{outdir}/{prefix}.tree", "-t", f"{fnb}-{tool}"])
          fi
        """)

        post_fh.write(f'''
            if os.path.exists("{outdir}/{prefix}.tree.svg"):
                shutil.copy("{outdir}/{prefix}.tree.svg", "{tool_rep_dir}/{fnb}-{tool}.tree.svg")
        ''')

        # need to be done once.
        if idx == 1:
            for rank in ("genus", "species", "strain"):
                if not os.path.exists(f"{p_repdir}/heatmap_TOOL-{tool}.{rank}.pdf"):
                    post_fh.write(f'''
                        merge_list_specTaxa_by_tool.pl {p_outdir}/*/{tool}/*.list -p {fnb} --top {heatmap_top} -l {rank} -otu {tmpdir}/{tool}.{rank}.otu.txt -output {tmpdir}/{tool}.{rank}.heatmap.matrix;
                        biom convert -i {tmpdir}/{tool}.{rank}.otu.txt -o {tool_rep_dir}/{fnb}-{tool}-{rank}.biom --to-hdf5 --table-type='OTU table' --process-obs-metadata taxonomy 2>/dev/null;
                        heatmap_distinctZ_noClust_zeroRowAllow.py --maxv 100 -s {heatmap_scale} --in {tmpdir}/{tool}.{rank}.heatmap.matrix --out {p_repdir}/heatmap_TOOL-{tool}.{rank}.pdf  2>/dev/null;
                    ''')

        post_fh.write(")&\n")
    
    post_fh.write("\nwait\n")

    fnb_rep_dir = f"{p_repdir}/{idx}_{fnb}"

    post_fh.write(f"""
        echo "==> Generating Radar Chart...";
        convert_list2radarChart.pl --level genus   --outdir {p_repdir} --outprefix radarchart_DATASET {fnb_rep_dir}/*/*.list.txt &
        convert_list2radarChart.pl --level species --outdir {p_repdir} --outprefix radarchart_DATASET {fnb_rep_dir}/*/*.list.txt &
        convert_list2radarChart.pl --level strain  --outdir {p_repdir} --outprefix radarchart_DATASET {fnb_rep_dir}/*/*.list.txt &

        echo "==> Generating matrix for heatmap by tools...";
        merge_list_specTaxa_by_dataset.pl {fnb_rep_dir}/*/*.list.txt --top {heatmap_top} -l genus   --otu {tmpdir}/{fnb}.genus.otu.txt --output {tmpdir}/{fnb}.genus.heatmap.matrix & 
        merge_list_specTaxa_by_dataset.pl {fnb_rep_dir}/*/*.list.txt --top {heatmap_top} -l species --otu {tmpdir}/{fnb}.species.otu.txt --output {tmpdir}/{fnb}.species.heatmap.matrix &
        merge_list_specTaxa_by_dataset.pl {fnb_rep_dir}/*/*.list.txt --top {heatmap_top} -l strain  --otu {tmpdir}/{fnb}.strain.otu.txt --output {tmpdir}/{fnb}.strain.heatmap.matrix &

        wait 

        biom convert -i {tmpdir}/{fnb}.genus.otu.txt -o {p_repdir}/{fnb}.genus.biom --to-hdf5 --table-type='OTU table' --process-obs-metadata taxonomy 2>/dev/null &
        biom convert -i {tmpdir}/{fnb}.species.otu.txt -o {p_repdir}/{fnb}.species.biom --to-hdf5 --table-type='OTU table' --process-obs-metadata taxonomy 2>/dev/null &
        biom convert -i {tmpdir}/{fnb}.strain.otu.txt -o {p_repdir}/{fnb}.strain.biom --to-hdf5 --table-type='OTU table' --process-obs-metadata taxonomy 2>/dev/null &

        wait

        echo "==> Generating heatmaps...";
        heatmap_distinctZ_noClust_zeroRowAllow.py --maxv 100 -s {heatmap_scale} --in {tmpdir}/{fnb.genus}.heatmap.matrix --out {p_repdir}/heatmap_DATASET-{fnb.genus}.pdf --title {fnb.genus} 2>/dev/null &
        heatmap_distinctZ_noClust_zeroRowAllow.py --maxv 100 -s {heatmap_scale} --in {tmpdir}/{fnb.species}.heatmap.matrix --out {p_repdir}/heatmap_DATASET-{fnb.species}.pdf --title {fnb.species} 2>/dev/null &
        heatmap_distinctZ_noClust_zeroRowAllow.py --maxv 100 -s {heatmap_scale} --in {tmpdir}/{fnb.strain}.heatmap.matrix --out {p_repdir}/heatmap_DATASET-{fnb.strain}.pdf --title {fnb.strain} 2>/dev/null &
    """)

    print(f"\nwait\n")
    print(f'echo "[END #{idx} {fnb}]"\n\n')


hl_flag = f"--highlight_list={hl_list}" if hl_list is not None and os.path.exists(hl_list) else ""
print(f"\necho \"==> Generating TOP{rep_top} contamination report...\"")
print(f"convert_list2report.pl --top {rep_top} --list {filelist_out} --setting {opt['setting']} --output {p_outdir} > {summary_out} &")
print(f"\necho \"==> Generating resource usage report...\"")
print(f"uge_helper -l {logfile} > {res_usage_out} &")
print(f"\necho \"==> Producing report XLSX file...\"")
print(f"generate_xlsx_report.pl {hl_flag} --list {filelist_out} --setting {opt['setting']} --output {p_outdir} &")
print("wait")

if tools["system"]["RUN_POST_PROCESS"]:
    print("\n[POST PROCESS] Generate report...\n")
    os.system(f"bash {post_script}")
    print("\n[POST PROCESS] Done.\n")

# clean up output directory
if not opt["debug"]:
    os.system(f"rm -rf {p_outdir}/script")
    os.system(f"rm -rf {p_outdir}/temp")
    os.system(f"rm -rf {p_outdir}/*_*")
    os.system(f"rm -rf {p_seqdir}")
    # os.system(f"rm -rf {p_logdir}") # Uncomment this line if you want to remove p_logdir as well

###############################################################

def prep_sequence(file_info, tools):
    # Flags
    FASTQ, FASTQSE, FASTQPE, FASTA, FASTA_EXTRACT, SPLITRIM_DIR = 0,0,0,0,0,0

    for tool in tools.keys():
        cmd = tools[tool]['COMMAND']
        FASTQ = 1 if '%FASTQ%' in cmd else 0
        FASTQSE = 1 if '%FASTQSE%' in cmd else 0
        FASTQPE = 1 if '%FASTQPE%' in cmd else 0
        FASTA = 1 if '%FASTA%' in cmd else 0
        FASTA_EXTRACT = 1 if '%FASTA_EXTRACT%' in cmd else 0
        SPLITRIM_DIR = 1 if '%SPLITRIM_DIR%' in cmd else 0

    FASTA = 1 if FASTA_EXTRACT else 0

    for count in sorted(file_info.keys()):
        print("\n[PREP_SEQ] Processing #{} sequence...".format(count))

        fnb = file_info[count]['PREFIX']
        if not fnb:
            if file_info[count].get('FASTQ'):
                fnb = file_info[count]['FASTQ'].split('/')[-1].split('.')[0]
            elif file_info[count].get('FASTQPE'):
                fnb = file_info[count]['FASTQPE'].split('/')[-1].split('.')[0]
            else:
                fnb = "Dataset{}".format(count)
            file_info[count]['PREFIX'] = fnb
            print("[PREP_SEQ] Generate filename base for output prefix: {}".format(fnb))
        else:
            print("[PREP_SEQ] PREFIX: {}".format(file_info[count]['PREFIX']))

        if not os.path.exists(file_info[count]['FASTQ']) and FASTQ:
            print("[PREP_SEQ] FASTQ not found! Generate FASTQ seq: {}/{}.fastq".format(p_seqdir, fnb))
            pe = file_info[count]['FASTQPE'].split(',')
            if os.path.exists(pe[0]) and os.path.exists(pe[1]):
                cmd = "cat {} {} > {}/{}.fastq".format(pe[0], pe[1], p_seqdir, fnb)
                os.system(cmd)
                file_info[count]['FASTQ'] = "{}/{}.fastq".format(p_seqdir, fnb)
        else:
            print("[PREP_SEQ] FASTQ: {}".format(file_info[count]['FASTQ']))

        file = file_info[count]['FASTQ']
        path = file.rsplit('/', 1)[0] if '/' in file else '.'


        # splitrim
        # TRIM_FIXL = 30
        # TRIM_MINQ = 20
        # if SPLITRIM_DIR:
        #     if not os.path.exists(f"{p_seqdir}/splitrim_fixL{TRIM_FIXL}Q{TRIM_MINQ}/{fnb}_splitrim.fastq"):
        #         _verbose(f"[PREP_SEQ] SPLITRIM_DIR not found! Generate SPLITRIM_DIR: {p_seqdir}/splitrim_fixL{TRIM_FIXL}Q{TRIM_MINQ}/")
        #         fastq = file_info[count]['FASTQ']
        #         cmd = f"{os.environ['EDGE_HOME']}/thirdParty/gottcha/bin/splitrim --inFile={fastq} --fixL={TRIM_FIXL} --recycle --minQ={TRIM_MINQ} --prefix={fnb} --outPath={p_seqdir}/splitrim_fixL{TRIM_FIXL}Q{TRIM_MINQ} > /dev/null 2>&1"
        #         exe = os.system(cmd)
        #         if exe:
        #             raise Exception("Can't not generate splitrim directory.")
        #         file_info[count]['SPLITRIM_DIR'] = f"{p_seqdir}/splitrim_fixL{TRIM_FIXL}Q{TRIM_MINQ}"
        #     else:
        #         _verbose(f"[PREP_SEQ] SPLITRIM_DIR: {p_seqdir}/splitrim_fixL{TRIM_FIXL}Q{TRIM_MINQ}")
        #     file_info[count]['SPLITRIM_DIR'] = f"{p_seqdir}/splitrim_fixL{TRIM_FIXL}Q{TRIM_MINQ}"


        # convert FASTQ to FASTA
        if 'FASTA' not in file_info[count] and FASTA:
            file_info[count]['FASTA'] = f"{path}/{fnb}.fa" if os.path.exists(f"{path}/{fnb}.fa") else None
            file_info[count]['FASTA'] = f"{path}/{fnb}.fasta" if os.path.exists(f"{path}/{fnb}.fasta") else None
            file_info[count]['FASTA'] = f"{p_seqdir}/{fnb}.fa" if os.path.exists(f"{p_seqdir}/{fnb}.fa") else None
            file_info[count]['FASTA'] = f"{p_seqdir}/{fnb}.fasta" if os.path.exists(f"{p_seqdir}/{fnb}.fasta") else None

            if file_info[count]['FASTA'] is None:
                _verbose(f"[PREP_SEQ] FASTA not found. Generate FASTA sequence: {p_seqdir}/{fnb}.fasta")
                # generate fasta
                cmd = f"fastq_to_fasta_fast < {file} > {p_seqdir}/{fnb}.fasta"
                os.system(cmd)
                file_info[count]['FASTA'] = f"{p_seqdir}/{fnb}.fasta"
            else:
                _verbose(f"[PREP_SEQ] FASTA found: {file_info[count]['FASTA']}")
        else:
            _verbose(f"[PREP_SEQ] FASTA: {file_info[count]['FASTA']}")


        # extract FASTA
        if 'FASTA_EXTRACT' not in file_info[count] and FASTA_EXTRACT:
            file_info[count]['FASTA_EXTRACT'] = f"{path}/{fnb}.extract.fa" if os.path.exists(f"{path}/{fnb}.extract.fa") else None
            file_info[count]['FASTA_EXTRACT'] = f"{path}/{fnb}.extract.fasta" if os.path.exists(f"{path}/{fnb}.extract.fasta") else None
            file_info[count]['FASTA_EXTRACT'] = f"{p_seqdir}/{fnb}.extract.fa" if os.path.exists(f"{p_seqdir}/{fnb}.extract.fa") else None
            file_info[count]['FASTA_EXTRACT'] = f"{p_seqdir}/{fnb}.extract.fasta" if os.path.exists(f"{p_seqdir}/{fnb}.extract.fasta") else None

            if file_info[count]['FASTA_EXTRACT'] is None:
                # generate fasta
                _verbose(f"[PREP_SEQ] Extracted FASTA not found. Extracting FASTA ({extract}): {p_seqdir}/{fnb}.extract.fasta")
                cmd = f"extract_random_sequences.pl -n {extract} -i {file_info[count]['FASTA']} -o {p_seqdir}/{fnb}.extract.fasta"
                _log(f"[PREP_SEQ] CMD={cmd}")
                os.system(cmd)
                file_info[count]['FASTA_EXTRACT'] = f"{p_seqdir}/{fnb}.extract.fasta"
            else:
                _verbose(f"[PREP_SEQ] FASTA_EXTRACT found: {file_info[count]['FASTA_EXTRACT']}")
        else:
            _verbose(f"[PREP_SEQ] FASTA_EXTRACT: {file_info[count]['FASTA_EXTRACT']}")

def param_replace(cmd, info, tool_ref, idx, tool):
    # variables defined in tool sections
    for key, val in tool_ref[tool].items():
        cmd = cmd.replace(f"%{key}%", val)

    # variables defined in [system]
    system = tool_ref["system"]
    for key, val in system.items():
        cmd = cmd.replace(f"%{key}%", val)

    # pre-defined variables
    cmd = cmd.replace("%PREFIX%", info[idx]["PREFIX"])
    cmd = cmd.replace("%FASTQ%", info[idx]["FASTQ"])
    cmd = cmd.replace("%FASTQSE%", info[idx]["FASTQSE"])
    cmd = cmd.replace("%FASTQPE%", info[idx]["FASTQPE"])
    cmd = cmd.replace("%FASTA%", info[idx]["FASTA"])
    cmd = cmd.replace("%SPLITRIM_DIR%", info[idx]["SPLITRIM_DIR"])
    cmd = cmd.replace("%FASTA_EXTRACT%", info[idx]["FASTA_EXTRACT"])
    cmd = cmd.replace("%TOOL%", tool)
    cmd = cmd.replace("%SERIAL%", str(idx))

    return cmd

def restore_settings(file):
    set = {}
    section = None
    count = 0
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            if line == '' or line.startswith('#') or line.startswith(';;'):
                continue
            if line.startswith('[') and line.endswith(']'):  # new section
                section = line[1:-1]
                if section in set:
                    raise ValueError(f'Section "{section}" already exists. Please check the setting file.')
                set[section] = {'ORDER': count}
                count += 1
                continue
            key, val = line.split('=', 1)
            key = key.strip().upper()
            set[section][key] = val.strip()
    return set

def _log(msg):
    if log_fh:
        log_fh.write(msg)

def _verbose(msg):
    if opt['verbose']:
        print(msg)
    _log(msg)

def _notify(msg):
    print(msg)
    _log(msg)

def _notifyError(msg):
    _log(msg)
    raise ValueError(msg)

def countFastq_exe(file):
    seq_count = 0
    total_length = 0
    with open(file, 'r') as fh:
        for line in fh:
            id = line
            seq = next(fh)
            seq = seq.strip()
            q_id = next(fh)
            q_seq = next(fh)
            length = len(seq)
            seq_count += 1
            total_length += length
    return (seq_count, total_length)

def getCpuUsage(pid=None):
    if pid is None:
        pid = str(os.getpid())
    pstree = subprocess.check_output(['pstree', '-lp', pid]).decode()
    pids = re.findall(r'\((\d+)\)', pstree)
    cpu = 0
    for pid in pids:
        usage = subprocess.check_output(['ps', '-p', pid, '-o', 'c']).decode().strip()
        if usage.isdigit():
            cpu += int(usage)
    return cpu / 100

def timeInterval(now):
    now = time.time() - now
    hours = int(now / 3600)
    minutes = int((now % 3600) / 60)
    seconds = int(now % 60)
    return f"{hours:02d}:{minutes:02d}:{seconds:02d}"

def usage():
    print("""

USAGE: python sample_contamination.py -s [INI_FILES] [INPUT FASTQ] ([INPUT FASTQ2] [INPUT FASTQ3]...) (-l [FILELIST]) 

    [INPUT FASTQ]        One or more input FASTAQ file need to be provided. Wildcard is allowed. User can also
	                     provide an filelist. Check \"-l\" option for more detail.

    -l [FILELIST]        This argument is optional. The filelist can be a simple FASTQ path list or in a certain
	                     format like /users/218817/bin/sample_contamination.filelist.txt.

    -s [SETTINGS_FILE]   A setting file with certain format is required for the wrapper. An example ini file can be 
                         found in /users/218817/bin/sample_contamination.settings.ini.
 
    -c [NUM]             Number of CPU (default: 4)
    
    -d                   Debug mode. Keep output and all temporary files.

    -o [OUTPUT_DIRECTORY]  

=================================================================================================================

Example:

python sample_contamination.py -s sample_contamination.settings.ini -o example_out sequence/*.fastq  
python sample_contamination.py -s sample_contamination.settings.ini -o example_out -l sequence_list.txt
python sample_contamination.py -s settings.txt -c 4 -o example_out seq1.fastq seq2.fastq seq3.fastq

""")
    sys.exit(1)