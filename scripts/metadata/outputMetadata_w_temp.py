import sys
import os
from os.path import dirname
from html_template import Template
from datetime import datetime

out_dir = sys.argv[1]
html_outfile = sys.argv[2]
projname = sys.argv[3]
out_dir_parts = out_dir.split('/')
projid = out_dir_parts[-1]

# Instantiate the variables
vars = {}

def output_html():
    vars["OUTPUTDIR"] = out_dir
    vars.setdefault("PROJNAME", projname)
    vars.setdefault("PROJID", projid)

    template = Template(filename=os.path.join(dirname(__file__), "edge_html_metadata.tmpl"), strict=False,
                        die_on_bad_params=False)
    template.set_params(vars)

    os.makedirs(os.path.join(out_dir, "Metadata"), exist_ok=True)

    with open(html_outfile, "w") as htmlfh:
        htmlfh.write(template.output())

def pull_sampleMetadata():
    metadata = f"{out_dir}/metadata_sample.txt"
    if os.path.exists(metadata):
        with open(metadata, 'r') as conf:
            for line in conf:
                line = line.strip()
                if not line.startswith('#'):
                    match = re.match(r'(.*)=(.*)', line)
                    if match:
                        key = match.group(1)
                        value = match.group(2)
                        if key == "study_title":
                            vars["SMD_STUDY_TITLE"] = value
                        elif key == "study_type":
                            vars["SMD_STUDY_TYPE"] = value
                        elif key == "sample_name":
                            vars["SMD_NAME"] = value
                        elif key == "sample_type":
                            vars["SMD_TYPE"] = value
                        elif key == "gender":
                            vars["OUT_SMD_GENDER"] = 1
                            vars["SMD_GENDER"] = value
                        elif key == "age":
                            vars["OUT_SMD_AGE"] = 1
                            vars["SMD_AGE"] = value
                        elif key == "host":
                            vars["OUT_SMD_HOST"] = 1
                            vars["SMD_HOST"] = value
                        elif key == "host_condition":
                            vars["OUT_SMD_HOST_CONDITION"] = 1
                            vars["SMD_HOST_CONDITION"] = value
                        elif key == "isolation_source":
                            vars["SMD_SOURCE"] = value
                        elif key == "collection_date":
                            vars["SMD_COLLECTION_DATE"] = value
                        elif key == "location":
                            vars["SMD_LOCATION"] = value
                        elif key == "city":
                            vars["SMD_CITY"] = value
                        elif key == "state":
                            vars["SMD_STATE"] = value
                        elif key == "country":
                            vars["SMD_COUNTRY"] = value
                        elif key == "lat":
                            vars["SMD_LAT"] = value
                        elif key == "lng":
                            vars["SMD_LNG"] = value
                        elif key == "experiment_title":
                            vars["SMD_EXP_TITLE"] = value
                        elif key == "sequencing_center":
                            vars["SMD_SEQ_CENTER"] = value
                        elif key == "sequencer":
                            vars["SMD_SEQUENCER"] = value
                        elif key == "sequencing_date":
                            vars["SMD_SEQ_DATE"] = value
        conf.close()

        if vars.get("SMD_TYPE") and vars["SMD_TYPE"] == "human":
            vars["SMD_HOST"] = ""

        if vars.get("SMD_TYPE"):
            type = vars["SMD_TYPE"]
            setSampleType(type)

        if vars.get("SMD_GENDER"):
            gender = vars["SMD_GENDER"]
            setGender(gender)
        else:
            setGender("male")

        if vars.get("SMD_HOST_CONDITION"):
            host_condition = vars["SMD_HOST_CONDITION"]
            setHostCondition(host_condition)
        else:
            setHostCondition("unknown")
    else:
        setSampleType("human")
        setGender("male")
        setHostCondition("unknown")
        
def set_sample_type(type):
    types = ['animal', 'environmental', 'human']
    for i in range(len(types)):
        item = {}
        item['SMD_TYPE'] = types[i]
        item['SMD_TYPE_LABEL'] = types[i].capitalize()
        item['SMD_TYPE_ID'] = f"metadata-sample-type{i}"
        if type == types[i]:
            item['SMD_TYPE_CHECKED'] = "checked"
        vars['LOOP_SMD_TYPES'].append(item)

def set_gender(gender):
    genders = ['male', 'female']
    for i in range(len(genders)):
        item = {}
        item['SMD_GENDER'] = genders[i]
        item['SMD_GENDER_LABEL'] = genders[i].capitalize()
        item['SMD_GENDER_ID'] = f"metadata-host-gender{i}"
        if gender == genders[i]:
            item['SMD_GENDER_CHECKED'] = "checked"
        vars['LOOP_SMD_GENDERS'].append(item)

def set_host_condition(type):
    types = ['healthy', 'diseased', 'unknown']
    for i in range(len(types)):
        item = {}
        item['SMD_HOST_CONDITION'] = types[i]
        item['SMD_HOST_CONDITION_LABEL'] = types[i].capitalize()
        item['SMD_HOST_CONDITION_ID'] = f"metadata-host-condition{i}"
        if type == types[i]:
            item['SMD_HOST_CONDITION_CHECKED'] = "checked"
        vars['LOOP_SMD_HOST_CONDITIONS'].append(item)

def pull_other_metadata():
    metadata = f"{out_dir}/metadata_other.txt"
    if os.path.exists(metadata):
        with open(metadata, 'r') as fh:
            vars['SMD_OTHER'] = fh.read()
