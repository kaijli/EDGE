#!/usr/bin/env python
import os
import argparse
import xlsxwriter
from os.path import dirname
from pathlib import Path

um = None  # user management
out = None
project_dir_names = None
usage = """
Usage: {0}
	Required:
		-out             output file
		-projects        project_dir_names, separated by comma
""".format(__file__)

parser = argparse.ArgumentParser(description='Convert Perl script to Python')
parser.add_argument('-um', type=str, help='user management')
parser.add_argument('-out', type=str, help='output file')
parser.add_argument('-projects', type=str, help='project_dir_names')
parser.add_argument('-help', action='store_true', help='display help message')

args = parser.parse_args()
um = args.um
out = args.out
project_dir_names = args.projects

if not project_dir_names and not out:
    print(usage)
    exit()

Path(dirname(out)).mkdir(parents=True, exist_ok=True)

# create sheets and write headers
n1, n2, n3, n4 = 1, 1, 1, 1
workbook = xlsxwriter.Workbook(out)
header_format = workbook.add_format({'bold': True})

if um:
    header1 = ["Project/Run Name", "Project/Run Owner", "Study Title", "Study Type", "Sample Name", "Sample Type",
               "Host", "Host Condition", "Gender", "Age", "Isolation Source", "Collection Date", "Location", "City",
               "State", "Country", "Lat", "Lng", "Experiment Title", "Sequencing Center", "Sequencer", "Sequencing Date"]
    header2 = ["Project/Run Name", "Project/Run Owner", "Date From", "Date To", "Location", "City", "State", "Country",
               "Lat", "Lng"]
    header3 = ["Project/Run Name", "Project/Run Owner", "Category", "Symptom"]
    header4 = ["Project/Run Name", "Project/Run Owner", "Field", "Value"]
else:
    header1 = ["Project/Run Name", "Study Title", "Study Type", "Sample Name", "Sample Type", "Host",
               "Host Condition", "Gender", "Age", "Isolation Source", "Collection Date", "Location", "City", "State",
               "Country", "Lat", "Lng", "Experiment Title", "Sequencing Center", "Sequencer", "Sequencing Date"]
    header2 = ["Project/Run Name", "Date From", "Date To", "Location", "City", "State", "Country", "Lat", "Lng"]
    header3 = ["Project/Run Name", "Category", "Symptom"]
    header4 = ["Project/Run Name", "Field", "Value"]

worksheet1 = workbook.add_worksheet("sample_metadata")
worksheet1.write_row("A{0}".format(n1), header1, header_format)
n1 += 1

worksheet2 = workbook.add_worksheet("travels")
worksheet2.write_row("A{0}".format(n2), header2, header_format)
n2 += 1

worksheet3 = workbook.add_worksheet("symptoms")
worksheet3.write_row("A{0}".format(n3), header3, header_format)
n3 += 1

worksheet4 = workbook.add_worksheet("other")
worksheet4.write_row("A{0}".format(n4), header4, header_format)
n4 += 1

# write metadata to sheets
for proj_dir in project_dir_names.split(","):
    confFile = f"{proj_dir}/config.txt"
    metadataFile = f"{proj_dir}/metadata_sample.txt"
    conf = getParams(confFile)
    metadata = getParams(metadataFile)
    travelsFile = f"{proj_dir}/metadata_travels.txt"
    symptomsFile = f"{proj_dir}/metadata_symptoms.txt"
    otherFile = f"{proj_dir}/metadata_other.txt"

    proj_name = conf['projname']
    owner = conf['projowner']
    if os.path.exists(metadataFile):
        writeSampleMetadata(owner, proj_name, metadata)
        writeTravelsMetadata(owner, proj_name, travelsFile)
        writeSymptomsMetadata(owner, proj_name, symptomsFile)
        writeOtherMetadata(owner, proj_name, otherFile)

workbook.close() # assuming $workbook is already defined

def writeSampleMetadata(owner, proj, meta):
    global n1, worksheet1
    if um:
        row = [proj, owner, meta['study_title'], meta['study_type'], meta['sample_name'], meta['sample_type'], meta['host'], meta['host_condition'], meta['gender'], meta['age'], meta['isolation_source'], meta['collection_date'], meta['location'], meta['city'], meta['state'], meta['country'], meta['lat'], meta['lng'], meta['experiment_title'], meta['sequencing_center'], meta['sequencer'], meta['sequencing_date']]
    else:
        row = [proj, meta['study_title'], meta['study_type'], meta['sample_name'], meta['sample_type'], meta['host'], meta['host_condition'], meta['gender'], meta['age'], meta['isolation_source'], meta['collection_date'], meta['location'], meta['city'], meta['state'], meta['country'], meta['lat'], meta['lng'], meta['experiment_title'], meta['sequencing_center'], meta['sequencer'], meta['sequencing_date']]
    
    worksheet1.write(f"A{n1}", row)
    n1 += 1

def writeTravelsMetadata(owner, proj, file):
    if not os.path.exists(file):
        return
    
    with open(file, "r") as f:
        lines = f.readlines()
        from_date = ""
        to_date = ""
        location = ""
        city = ""
        state = ""
        country = ""
        lat = ""
        lng = ""
        for line in lines:
            line = line.strip()
            if line.startswith("#"):
                continue
            parts = line.split("=")
            if len(parts) == 2:
                key = parts[0].strip()
                value = parts[1].strip()
                if key == "travel-date-from":
                    from_date = value
                elif key == "travel-date-to":
                    to_date = value
                elif key == "travel-location":
                    location = value
                elif key == "city":
                    city = value
                elif key == "state":
                    state = value
                elif key == "country":
                    country = value
                elif key == "lat":
                    lat = value
                elif key == "lng":
                    lng = value
                    row = []
                    if um:
                        row = [proj, owner, from_date, to_date, location, city, state, country, lat, lng]
                    else:
                        row = [proj, from_date, to_date, location, city, state, country, lat, lng]
                    worksheet2.write("A" + str(n2), row)
                    n2 += 1

def writeSymptomsMetadata(owner, proj, file):
    if not os.path.exists(file):
        return
    
    with open(file, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) == 2:
                key = parts[0].strip()
                value = parts[1].strip()
                row = []
                if um:
                    row = [proj, owner, key, value]
                else:
                    row = [proj, key, value]
                worksheet3.write("A" + str(n3), row)
                n3 += 1

def writeOtherMetadata(owner, proj, file, um, worksheet4, n4):
    if not os.path.exists(file):
        return

    with open(file, 'r') as sm:
        for line in sm:
            line = line.strip()
            if line.startswith('#'):
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                row = []
                if um:
                    row = [proj, owner, key, value]
                else:
                    row = [proj, key, value]
                worksheet4.write(f"A{n4}", row)
                n4 += 1


def getParams(config):
    sys = {}
    if os.path.exists(config):
        with open(config, 'r') as conf:
            for line in conf:
                line = line.strip()
                if line.startswith('#'):
                    continue
                if '=' in line:
                    key, value = line.split('=', 1)
                    sys[key] = value
    return sys
