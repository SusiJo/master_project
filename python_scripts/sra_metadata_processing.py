#!/usr/bin/python3

# This script processes metadata xml, csv files downloaded from SRA

# imports
import click
import os.path
import sys
import logging
import time
from pathlib import Path
import pandas as pd
import re
import xml.etree.ElementTree as ET
import fileUtils.file_handling as fh

# Create logger
logger = logging.getLogger('SRA metadata reader')
# Create console handler
ch = logging.StreamHandler()
# Create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)
logger.setLevel(logging.INFO)


@click.command()
@click.option('-x', '--xml', prompt='path to xml metadata files',
              help='Path to folder with metadata in xml format',
              required=True)
@click.option('-c', '--csv', prompt='path to csv metadata files',
              help='Path to folder with metadata in csv format',
              required=True)
@click.option('--csv_only/--no_csv_only', prompt='write all csv to table?',
              help='If csv_only specified, additionally writes all csv files to one output table',
              required=False)
@click.option('-o', '--outpath', prompt='path for output table (Experiment & Runinfo)',
              help='Table with metadata (Experiment xml, Runinfo csv) metadata in csv format',
              required=True)
def main(xml, csv, outpath, csv_only):
    start_time = time.time()
    # printed to STDOUT on command line
    logger.info('Get files')
    all_xml_files = fh.get_files(xml, '*.xml')
    all_csv_files = fh.get_files(csv, '*.csv')

    logger.info('Parsing xml files...')
    metadata_xml = parse_xml(all_xml_files, all_csv_files)

    logger.info('Writing output to table...')
    write_xml_table(outpath, metadata_xml)

    if csv_only:
        logger.info('Writing all csv files into table...')
        parse_all_fields_csv(all_csv_files, outpath)

    end_time = time.time()
    logger.info('Process finished in ' + str(round(end_time - start_time, 2)) + "sec")


def parse_all_fields_csv(files, outpath):
    """Reads all csv files and produces one table with all fields

    :param files
    :param outpath
    """
    p = Path(outpath)
    ext = p.suffix
    new_outpath = p.rename(Path(p.parent, "metadata_all_csv" + ext))
    dfs = list()
    for filename in files:
        df = pd.read_csv(filename)
        dfs.append(df)

    frame = pd.concat(dfs, axis=0, ignore_index=True)
    frame.to_csv(new_outpath, header=True, index=False)


def parse_xml(xml_files, csv_files):
    """Reads xml files which are not consistent

    :param xml_files: as list
    :param csv_files: as list
    :return: metadata: dict
    """
    metadata = {}

    for i, file in enumerate(xml_files):
        tree = ET.parse(file)
        root = tree.getroot()
        srr = os.path.basename(file).split(".")[0]
        source_name = ''
        body_site = ''
        subject_status = ''
        project = ''
        biosample = ''
        gender = ''
        study_disease = ''
        histological_type = ''
        is_tumor = ''
        age = ''
        tissue = ''
        phenotype = ''
        library_strat = ''
        organism = ''
        # print(srr)
        for h, child in enumerate(root):
            for j, ch in enumerate(child):
                for idx, c in enumerate(ch):
                    for k, e in enumerate(c):
                        if e.tag == 'SCIENTIFIC_NAME':
                            organism = root[h][j][idx][k].text
                            # print(organism)
                        for i, s in enumerate(e):
                            # print(s.attrib)
                            if s.tag == 'LIBRARY_STRATEGY':
                                library_strat = root[h][j][idx][k][i].text
                                # print(library_strat)
                            if s.tag == 'EXTERNAL_ID' and s.get('namespace') == 'BioProject':
                                # print(s.attrib) #{'namespace': 'BioProject'}
                                project = root[h][j][idx][k][i].text
                                # print(project)
                            if s.tag == 'EXTERNAL_ID' and s.get('namespace') == 'BioSample':
                                biosample = root[h][j][idx][k][i].text
                                # print(biosample)
                            if s.tag == 'TAG' and s.text == 'source_name':       # body_site
                                source_name = root[0][j][idx][k][i + 1].text
                                # print(body_site)
                            if s.tag == 'TAG' and s.text == 'subject status':
                                subject_status = root[0][j][idx][k][i + 1].text
                                # print(subject_status)
                            if s.tag == 'TAG' and (s.text == 'gender' or s.text == 'sex'):
                                gender = root[0][j][idx][k][i + 1].text
                                # print(gender)
                            if s.tag == 'TAG' and s.text == 'study disease':
                                study_disease = root[0][j][idx][k][i + 1].text
                            if s.tag == 'TAG' and s.text == 'histological type':
                                histological_type = root[0][j][idx][k][i + 1].text
                            if s.tag == 'TAG' and s.text == 'is tumor':
                                is_tumor = root[0][j][idx][k][i + 1].text
                                # print(is_tumor)
                            if s.tag == 'TAG' and s.text == 'age':
                                age = root[0][j][idx][k][i + 1].text
                                # print(age)
                            if s.tag == 'TAG' and s.text == 'tissue':
                                tissue = root[0][j][idx][k][i + 1].text
                                # print(tissue)
                            if s.tag == 'TAG' and s.text == 'phenotype':
                                phenotype = root[0][j][idx][k][i + 1].text
                                # print(phenotype)
                            if s.tag == 'TAG' and s.text == 'body site':
                                body_site = root[0][j][idx][k][i + 1].text
                                # print(body_site)
                            if source_name == '':
                                source_name = body_site
                            # if project not at that place, keep searching
                            if project == '':
                                if s.tag == 'TAG' and s.text == 'parent_bioproject':
                                    project = root[0][j][idx][k][i + 1].text
                            # if still empty
                            if project == '':
                                for o, m in enumerate(s):
                                    if m.tag == 'LABEL':
                                        if re.search("PRJ", m.text):
                                            project = m.text
                if project == '' and ch.tag == 'STUDY':
                    if re.search("PRJ", ch.get('alias')):
                        # print(ch.get('alias'))
                        project = ch.get('alias')

        # fill metadata dict
        if srr not in metadata:
            metadata[srr] = [project, biosample, library_strat, organism, gender, age, source_name, subject_status,
                             tissue, phenotype, is_tumor, study_disease, histological_type]

    # look-up condition in csv file for appending this info
    csv_dict = fh.parse_csv(csv_files)
    for key in metadata.keys():
        for k in csv_dict.keys():
            if key == k:
                # condition info from csv file
                metadata[key].append(csv_dict[k][1])
                # if project still empty, take from csv dict
                if metadata[key][0] == '':
                    metadata[key][0] = csv_dict[k][0]
    return metadata


def write_xml_table(outpath, metadata):
    """Write table with csv and xml info to table

    :param outpath: path to output table .csv
    :param metadata: dict
    """
    df = pd.DataFrame.from_dict(metadata, orient='index',
                                columns=['Project', 'BioSample', 'Library_Strategy', 'Organism', 'Gender', 'Age',
                                         'Source_name', 'Subject_status', 'Tissue', 'Phenotype', 'Is_Tumor',
                                         'StudyDisease', 'HistologicalType', 'Condition'])

    df.reset_index(level=0, inplace=True)
    df.rename({'index': 'RunID'}, axis=1, inplace=True)
    df.to_csv(outpath, header=True, index=False)


if __name__ == "__main__":
    sys.exit(main())
