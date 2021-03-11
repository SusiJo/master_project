# This script processes clinical metadata from ICGC from FILE and DONOR ENDPOINT in JSON

# imports
import click
import glob
import os.path
import sys
import logging
import time
import json
from dictor import dictor
import pandas as pd

# Create logger
logger = logging.getLogger('Json parser')
# Create console handler
ch = logging.StreamHandler()
# Create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)
logger.setLevel(logging.INFO)


@click.command()
@click.option('-f', '--file_endpt', prompt='path to folder with ICGC FILE metadata',
              help='Path to folder with ICGC FILE metadata',
              required=True)
@click.option('-d', '--donor_endpt', prompt='path to folder with ICGC DONOR metadata',
              help='Path to folder with ICGC DONOR metadata',
              required=True)
@click.option('-o', '--outpath', prompt='output table',
               help='Table with metadata in csv format',
               required=True)
def main(file_endpt, donor_endpt, outpath):
    start_time = time.time()

    # printed to STDOUT on command line
    logger.info('Get json files FILE endpoint')
    allfiles = get_files(file_endpt, '*.json')

    logger.info('Get json files DONOR endpoint')
    donorfiles = get_files(donor_endpt, '*.json')

    logger.info('Parse json files')
    metadata_allfiles = parse_icgc_json_files(allfiles, donorfiles)

    logger.info('Write metadata to table')
    table = write_icgc_table(outpath, metadata_allfiles)

    end_time = time.time()
    logger.info('Process finished in ' + str(round(end_time - start_time, 2)) + "sec")


def get_files(inpath, ext):
    """Create list of json files
    @:param: path to files
    @:param: file extension, json
    @:return: list of files
    """
    os.chdir(inpath)
    files = [os.path.abspath(os.path.basename(f)) for f in glob.glob(inpath + ext)]
    return files


def parse_icgc_json_files(files, donors):
    """Read json files from FILE & DONOR ENDPOINT and extract relevant values
    @:param: list of paths to files, FILE ENDPOINT
    @:param: list of paths  to files, DONOR ENDPOINT
    @:return: metadata dictionary
    """

    metadata = {}
    specimen_ids = []

    """Parse FILE ENDPOINT info store in metadata dict"""
    for f in files:
        # split absolute path to get basename
        file = f.split('/')[-1]
        with open(f) as json_file:
            data = json.load(json_file)
            # print(data)
            file_id = os.path.basename(file).split(".")[0]
            file_name = os.path.basename(file).split(".")[1]
            donor_id = dictor(data, "donors.0.donorId")
            specimen_type = dictor(data, "donors.0.specimenType.0")
            # print(dictor(data, "donors.0.specimenType.0"))
            sample_id = dictor(data, "donors.0.sampleId.0")
            specimen_id = dictor(data, "donors.0.specimenId.0")
            # print(file_id, file_name, donor_id)
            specimen_ids.append(specimen_id)

            if file_id not in metadata.keys():
                metadata[file_id] = [file_name, donor_id, sample_id, specimen_id, specimen_type]

    """Parse DONOR ENDPOINT info and append to metadata dict"""
    for d in donors:
        # print(d)
        donor = d.split('/')[-1]
        with open(donor) as json_file:
            data = json.load(json_file)
            sample_type = ''

            file_id = os.path.basename(donor).split(".")[0]
            # file_name = os.path.basename(donor).split(".")[1]
            # donor_id = os.path.basename(donor).split(".")[2]
            gender = dictor(data, "gender")
            vital_status = dictor(data, "vitalStatus")
            age_at_index = dictor(data, "ageAtDiagnosis")
            specimen = dictor(data, "specimen")
            tumor_subtype = dictor(data, "tumourSubtype")
            primary_site = dictor(data, "primarySite")
            primary_diagnosis = dictor(data, "tumourType")
            survival_time = dictor(data, "survivalTime")
            tumor_stage = dictor(data, "tumourStageAtDiagnosis")
            diagnosisIcd10 = dictor(data, "diagnosisIcd10")
            # match file and donor specimen
            for s in specimen_ids:
                for i, e in enumerate(specimen):
                    for k, v in e.items():
                        # print(s, i, e, k, v)
                        if v == s:
                            # print(e['type'])
                            # extract sample type
                            sample_type = e['type']
                            # primary_diagnosis = e['tumourHistologicalType'] # not found for some!!

            for k, v in metadata.items():
                if k == file_id:
                    v.append(sample_type)
                    v.append(primary_site)
                    v.append(primary_diagnosis)
                    v.append(tumor_subtype)
                    v.append(gender)
                    v.append(vital_status)
                    v.append(age_at_index)
                    v.append(survival_time)
                    v.append(tumor_stage)
                    v.append(diagnosisIcd10)

    return metadata


def write_icgc_table(outpath, metadata):
    """Write metadata to table
    @:param: path to store output table, CSV FORMAT
    @:param: metadata dictionary
    @:return: metadata table
    """
    table = pd.DataFrame.from_dict(metadata, orient='index',
                                   columns=['File_Name', 'Donor_ID', 'Sample_Id', 'Specimen_Id', 'Specimen_Type',
                                            'Sample_Type', 'Primary_Site', 'Primary_diagnosis', 'Tumor_Subtype', 'Gender',
                                            'Vital_status', 'Age_at_index', 'Survival_Time', 'Tumor_stage', 'Icd10'])

    table.reset_index(level=0, inplace=True)
    table.rename({'index': 'File_ID'}, axis=1, inplace=True)
    table.to_csv(outpath, index=False)
    return table


if __name__ == "__main__":
    sys.exit(main())
