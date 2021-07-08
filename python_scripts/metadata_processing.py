""" This script serves to process metadata from TCGA, ICGC, SRA and outputs a standardized table
    Metadata Table [sample_id, case_id, condition, bio_project]

    usage:
    python metadata_processing.py -i <path-to-fileEndpt-json> ... -s <path-to-csv-files> -t <path-to-json>
    -t <path-to-json> ... -o <outpath.csv>
"""

# imports
import os
import click
import os.path
import sys
import re
import json
from dictor import dictor
import pandas as pd
import fileUtils.file_handling as fh


@click.command()
@click.option('-i', '--icgc', prompt='paths to ICGC json metadata files FILE Endpoint',
              help='Path to folder with ICGC metadata in json format', multiple=True)
@click.option('-s', '--sra', prompt='path to SRA csv metadata files',
              help='Path to folder with SRA metadata in csv format')
@click.option('-t', '--tcga', prompt='paths to TCGA json metadata files',
              help='Path to folder with TCGA metadata in json format', multiple=True)
@click.option('-o', '--outpath', prompt='output path for storing metadata table',
              help='Table with metadata in csv format [Sample_ID, Case_ID, Condition, Project]',
              required=True)
def main(icgc, sra, tcga, outpath):
    dfs = list()

    # ICGC (several folders can be given)
    for folder in icgc:
        # print(folder)
        icgc_files = fh.get_files(folder, '*.json')
        icgc_metadata = parse_icgc_json_files(icgc_files)
        df = pd.DataFrame.from_dict(icgc_metadata, orient='index',
                                 columns=['CaseID', 'SampleType', 'Project'])

        df.reset_index(level=0, inplace=True)
        df.rename({'index': 'FileID'}, axis=1, inplace=True)
        dfs.append(df)

    # SRA (assumed to be in one folder)
    sra_files = fh.get_files(sra, '*.csv')
    sra_metadata = fh.parse_csv(sra_files)
    df_sra = pd.DataFrame.from_dict(sra_metadata, orient='index',
                                columns=['CaseID', 'SampleType', 'Project'])

    df_sra.reset_index(level=0, inplace=True)
    df_sra.rename({'index': 'File_ID'}, axis=1, inplace=True)
    dfs.append(df_sra)

    # TCGA (several folders can be given)
    for folder in tcga:
        # print(folder)
        tcga_files = fh.get_files(folder, '*.json')
        tcga_metadata = parse_tcga_json_files(tcga_files)
        df = pd.DataFrame.from_dict(tcga_metadata, orient='index',
                                    columns=['CaseID', 'SampleType', 'Project'])

        df.reset_index(level=0, inplace=True)
        df.rename({'index': 'FileID'}, axis=1, inplace=True)
        dfs.append(df)

    # Make Dataframe
    frame = pd.concat(dfs, axis=0, ignore_index=True)
    frame.to_csv(outpath, index=False)


def parse_tcga_json_files(files):
    """Read json files from list and extract relevant values

    :param files: as list
    :return: tcga_metadata: dict
    """

    tcga_metadata = {}
    normal_pattern = re.compile('normal', re.IGNORECASE)

    for i, file in enumerate(files):
        file = file.split('/')[-1]
        with open(file) as json_file:
            data = json.load(json_file)
            file_id = os.path.basename(file).split(".")[0]
            case_id = dictor(data, "data.cases.0.case_id")
            submitter_id = dictor(data, "data.submitter_id")  # = file_name
            sample_type = dictor(data, "data.cases.0.samples.0.sample_type")
            project = dictor(data, "data.cases.0.project.project_id")
            if bool(re.search(normal_pattern, sample_type)):
                sample_type = 'normal'
            else:
                sample_type = 'tumor'

            if submitter_id not in tcga_metadata.keys():
                tcga_metadata[submitter_id] = [case_id, sample_type, project]
    return tcga_metadata


def parse_icgc_json_files(files):
    """Read json files from list and extract relevant values

    :param files: as list
    :return: icgc_metadata: dict
    """

    icgc_metadata = {}
    normal_pattern = re.compile('normal', re.IGNORECASE)

    for file in files:
        print(file)
        file = file.split('/')[-1]
        with open(file) as json_file:
            data = json.load(json_file)
            file_id = os.path.basename(file).split(".")[0]
            file_name = os.path.basename(file).split(".")[1]
            donor_id = dictor(data, "donors.0.donorId")
            specimen_type = dictor(data, "donors.0.specimenType.0")
            project = dictor(data, "donors.0.projectCode")
            if bool(re.search(normal_pattern, specimen_type)):
                specimen_type = 'normal'
            else:
                specimen_type = 'tumor'

            if file_name not in icgc_metadata.keys():
                icgc_metadata[file_name] = [donor_id, specimen_type, project]

    return icgc_metadata


if __name__ == "__main__":
    sys.exit(main())


