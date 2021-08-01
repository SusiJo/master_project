#!/usr/bin/python3

# This script processes metadata json files downloaded from TCGA
# wget 'https://api.gdc.cancer.gov/files/UUID?expand=cases,cases.demographic,cases.samples,annotations,cases.diagnoses,cases.diagnoses.treatments,cases.project,cases.project.program,analysis.metadata.read_groups&pretty=true'


# imports
import click
import re
import glob
import os.path
import sys
import re
import logging
import time
import json
from dictor import dictor
import pandas as pd
import fileUtils.file_handling as fh


# Create logger
logger = logging.getLogger('Metadata table creator')
# Create console handler
ch = logging.StreamHandler()
# Create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)
logger.setLevel(logging.INFO)


@click.command()
@click.option('-i', '--inpath', prompt='path to json metadata files',
              help='Path to folder with TCGA metadata in json format',
              required=True)
@click.option('-o', '--outpath', prompt='path to output table',
              help='Table with metadata in csv format',
              required=True)
def main(inpath, outpath):
    start_time = time.time()
    # printed to STDOUT on command line
    logger.info('Get json files')
    allfiles = fh.get_files(inpath, '*.json')

    logger.info('Parse json files')
    metadata_allfiles = parse_tcga_json_files(allfiles)

    logger.info('Write metadata to table')
    write_tcga_table(outpath, metadata_allfiles)

    end_time = time.time()
    logger.info('Process finished in ' + str(round(end_time - start_time, 2)) + "sec")


def parse_tcga_json_files(files):
    """Read json files from list and extract relevant values

    :param files: list of files
    :return: metadata dictionary
    """

    metadata = {}
    # normal_pattern = re.compile('normal', re.IGNORECASE)
    dead = re.compile('dead', re.IGNORECASE)
    tumor_stages = {'i': 1, 'ii': 2, 'iii': 3, 'iv': 4}

    for i, file in enumerate(files):
        file = file.split('/')[-1]
        with open(file) as json_file:
            data = json.load(json_file)
            tumor_stage = ''
            file_id = os.path.basename(file).split(".")[0]
            case_id = dictor(data, "data.cases.0.case_id")
            submitter_id = dictor(data, "data.submitter_id")  # = file_name
            # file_name = dictor(data, "data.file_name")
            sample_type = dictor(data, "data.cases.0.samples.0.sample_type")
            project = dictor(data, "data.cases.0.project.project_id")
            # print(sample_type)
            # if re.search(normal_pattern, sample_type):
            #     sample_type = "normal"
            # else:
            #     sample_type = "tumor"
            primary_diagnosis = dictor(data, "data.cases.0.diagnoses.0.primary_diagnosis")
            gender = dictor(data, "data.cases.0.demographic.gender")
            vital_status = dictor(data, "data.cases.0.demographic.vital_status")
            age_at_index = dictor(data, "data.cases.0.demographic.age_at_index")
            stage = dictor(data, "data.cases.0.diagnoses.0.tumor_stage").split(' ')[-1]
            # split stage and transform from 'iia' format to 2
            stage = re.split('[abc]', stage)[0]
            if stage in tumor_stages.keys():
                tumor_stage = tumor_stages.get(stage)

            primary_site = dictor(data, "data.cases.0.diagnoses.0.site_of_resection_or_biopsy")
            diagnosisIcd10 = dictor(data, "data.cases.0.diagnoses.0.icd_10_code")
            if re.search(dead, vital_status):
                survival_time = dictor(data, "data.cases.0.demographic.days_to_death")
            else:
                survival_time = ''

            # Save to dictionary
            if file_id not in metadata.keys():
                metadata[file_id] = [submitter_id, case_id, project, primary_site, sample_type, primary_diagnosis,
                                     gender, vital_status, age_at_index, survival_time, tumor_stage, diagnosisIcd10]

    return metadata


def write_tcga_table(outpath, metadata):
    """Write metadata to table

    :param outpath: to metadata_table .csv
    :param metadata: dict
    """
    table = pd.DataFrame.from_dict(metadata, orient='index',
                                   columns=['Submitter_ID', 'Case_ID', 'Project', 'Primary_site', 'Sample_type',
                                            'Primary_diagnosis', 'Gender', 'Vital_status', 'Age_at_index',
                                            'Survival_time', 'Tumor_stage', 'Icd10'])
    table.reset_index(level=0, inplace=True)
    table.rename({'index': 'File_ID'}, axis=1, inplace=True)
    table.to_csv(outpath, index=False)


if __name__ == "__main__":
    sys.exit(main())
