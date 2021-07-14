import os
import click
import logging
from bs4 import BeautifulSoup
import fileUtils.file_handling as fh
import pandas as pd
import re
import time
import sys

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
@click.option('-o', '--outpath', prompt='path for output table ',
              help='Table with rich metadata and short metadata table in csv format.\n'
                   'Extension is filled automatically.',
              required=True)
def main(xml, outpath):
    start_time = time.time()

    logger.info('Get files')
    all_xml_files = fh.get_files(xml, '*.xml')

    logger.info('Parsing xml files...')
    parse_soup(all_xml_files, outpath)

    end_time = time.time()
    logger.info('Process finished in ' + str(round(end_time - start_time, 2)) + "sec")


def parse_soup(xml_files, outpath):
    """
    Function to parse a folder with xml files from NCBI
    :param xml_files: path to folder
    :param outpath: path for metadata.csv sheet
    """

    sample_attributes_dict = {}
    metadata = {}
    sample_attributes_list = []

    for i, file in enumerate(xml_files):
        # define &  clear all fields for metadata sheet
        # experiment_title = ''
        # study_ref = ''
        library_name = ''
        # library_strategy = ''
        # library_selection = ''
        bio_sample = ''
        sample_title = ''
        taxon_id = ''
        taxon_name = ''
        age = ''
        gender = ''  # sex, Sex
        disease = ''  # study_disease
        body_site = ''  # body site, body_site
        source_name = ''
        histological_type = ''
        tissue_type = ''
        is_tumor = ''
        # disease_state = ''  # disease status
        subject_status = ''  # subject is affected
        hbv_infection = ''
        cirrhosis = ''  # primary tumor cirrhosis
        tumor_stage = ''  # staging,  stage
        tumor_grade = ''  # total grade, primary tumor grade,
        cell_type = ''
        primary_tumor_cirrhosis = ''
        sample_type = ''
        clinical_characteristics = ''
        phenotype = ''
        tumor_type = ''
        tumor_status = ''
        cancer_type = ''
        diagnosis = ''
        health_state = ''
        label = ''
        description = ''
        rna_source = ''

        with open(file, 'r') as f:
            srr = os.path.basename(file).split(".")[0]
            soup = BeautifulSoup(f, "lxml-xml")
            # print(srr)
            experiment_title = soup.EXPERIMENT.TITLE.text

            study_ref = soup.EXPERIMENT.STUDY_REF.IDENTIFIERS.PRIMARY_ID.text

            if soup.DESIGN.LIBRARY_DESCRIPTOR.LIBRARY_NAME is not None:
                library_name = soup.DESIGN.LIBRARY_DESCRIPTOR.LIBRARY_NAME.text

            library_strategy = soup.DESIGN.LIBRARY_DESCRIPTOR.LIBRARY_STRATEGY.text

            library_selection = soup.DESIGN.LIBRARY_DESCRIPTOR.LIBRARY_SELECTION.text

            # if soup.PLATFORM.ILLUMINA.INSTRUMENT_MODEL is None:
            #     if soup.PLATFORM.ABI_SOLID.INSTRUMENT_MODEL is not None:
            #         instrument_model = soup.PLATFORM.ABI_SOLID.INSTRUMENT_MODEL.text

            # if soup.PLATFORM.ILLUMINA.INSTRUMENT_MODEL is not None:
            #     instrument_model = soup.PLATFORM.ILLUMINA.INSTRUMENT_MODEL.text
            #     print("instrument model:\t ", soup.PLATFORM.ILLUMINA.INSTRUMENT_MODEL.text)

            if soup.STUDY.IDENTIFIERS.EXTERNAL_ID is not None:
                bio_sample = soup.STUDY.IDENTIFIERS.EXTERNAL_ID.text
            if soup.SAMPLE.IDENTIFIERS.EXTERNAL_ID is not None:
                bio_sample = soup.SAMPLE.IDENTIFIERS.EXTERNAL_ID.text

            if soup.SAMPLE.TITLE is not None:
                sample_title = soup.SAMPLE.TITLE.text

            if soup.SAMPLE.SAMPLE_NAME.TAXON_ID is not None:
                taxon_id = soup.SAMPLE.SAMPLE_NAME.TAXON_ID.text

            if soup.SAMPLE.SAMPLE_NAME.SCIENTIFIC_NAME is not None:
                taxon_name = soup.SAMPLE.SAMPLE_NAME.SCIENTIFIC_NAME.text

            if soup.SAMPLE.SAMPLE_NAME.COMMON_NAME is not None:
                taxon_name = soup.SAMPLE.SAMPLE_NAME.COMMON_NAME.text

            if soup.SAMPLE.SAMPLE_ATTRIBUTES is not None:
                if soup.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.TAG.text == 'source_name':
                    source_name = soup.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.VALUE.text

                    # Extract all possible sample attribute tags
                sample_attributes = [(repr(sibling.TAG.text), repr(sibling.VALUE.text)) for sibling in
                                     soup.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.next_siblings]

                # Iterate over all sample attributes
                for sibling in soup.SAMPLE.SAMPLE_ATTRIBUTES.SAMPLE_ATTRIBUTE.next_siblings:
                    if sibling.TAG.text == 'age' or sibling.TAG.text == 'DONOR_AGE' or \
                            sibling.TAG.text == 'patient age at resection':
                        age = str("'" + sibling.VALUE.text)
                    if sibling.TAG.text == 'gender' or sibling.TAG.text == 'sex' \
                            or sibling.TAG.text == 'Sex' or sibling.TAG.text == 'patient sex' \
                            or sibling.TAG.text == 'DONOR_SEX':
                        gender = sibling.VALUE.text
                    if sibling.TAG.text == 'study disease' or sibling.TAG.text == 'disease':
                        disease = sibling.VALUE.text
                    if sibling.TAG.text == 'body site' or sibling.TAG.text == 'body_site' \
                            or sibling.TAG.text == 'organism_part' or sibling.TAG.text == 'OrganismPart' \
                            or sibling.TAG.text == 'source_name':
                        body_site = sibling.VALUE.text
                    if sibling.TAG.text == 'histological type':
                        histological_type = sibling.VALUE.text
                    if sibling.TAG.text == 'tissue type' or sibling.TAG.text == 'tissue type/source' \
                            or sibling.TAG.text == 'tissue' or sibling.TAG.text == 'tissue source' \
                            or sibling.TAG.text == 'tissue_type' or sibling.TAG.text == 'TISSUE_TYPE':
                        tissue_type = sibling.VALUE.text
                    if sibling.TAG.text == 'is tumor':
                        is_tumor = sibling.VALUE.text
                    if sibling.TAG.text == 'disease state' or sibling.TAG.text == 'disease status' \
                            or sibling.TAG.text == 'DISEASE' or sibling.TAG.text == 'infection status':
                        disease_state = sibling.VALUE.text
                    if sibling.TAG.text == 'subject status' or sibling.TAG.text == 'status':
                        subject_status = sibling.VALUE.text
                    if sibling.TAG.text == 'hbv infection':
                        hbv_infection = sibling.VALUE.text
                    if sibling.TAG.text == 'cirrhosis':
                        cirrhosis = sibling.VALUE.text
                    if sibling.TAG.text == 'tumor stage' or sibling.TAG.text == 'staging' \
                            or sibling.TAG.text == 'disease_stage' or sibling.TAG.text == 'stage' \
                            or sibling.TAG.text == 'disease stage':
                        tumor_stage = sibling.VALUE.text
                    if sibling.TAG.text == 'total grade' or sibling.TAG.text == 'tumor grade' \
                            or sibling.TAG.text == 'primary tumor grade':
                        tumor_grade = sibling.VALUE.text
                    if sibling.TAG.text == 'cell type':
                        cell_type = sibling.VALUE.text
                    if sibling.TAG.text == 'primary tumor cirrhosis':
                        primary_tumor_cirrhosis = sibling.VALUE.text
                    if sibling.TAG.text == 'sample_type' or sibling.TAG.text == "BIOMATERIAL_TYPE":
                        sample_type = sibling.VALUE.text
                    if sibling.TAG.text == 'clinical characteristics':
                        clinical_characteristics = sibling.VALUE.text
                    if sibling.TAG.text == 'phenotype' or sibling.TAG.text == 'Phenotype':
                        phenotype = sibling.VALUE.text
                    if sibling.TAG.text == 'tumor type':
                        tumor_type = sibling.VALUE.text
                    if sibling.TAG.text == 'tumor status':
                        tumor_status = sibling.VALUE.text
                    if sibling.TAG.text == 'cancer type':
                        cancer_type = sibling.VALUE.text
                    if sibling.TAG.text == 'diagnosis':
                        diagnosis = sibling.VALUE.text
                    if sibling.TAG.text == 'DONOR_HEALTH_STATUS' or sibling.TAG.text == 'health_state' \
                            or sibling.TAG.text == 'health state':
                        health_state = sibling.VALUE.text
                    if sibling.TAG.text == 'label':
                        label = sibling.VALUE.text
                    if sibling.TAG.text == 'Description':
                        description = sibling.VALUE.text
                    if sibling.TAG.text == 'rna source':
                        rna_source = sibling.VALUE.text

                    # save list of complete sample attribute possibilities for filtering interesing keywords
                    if repr(sibling.TAG.text) not in sample_attributes_dict.keys():
                        sample_attributes_dict[repr(sibling.TAG.text)] = i

            # save each Run accession and associated metadata
            if srr not in metadata.keys():
                metadata[srr] = [experiment_title, study_ref, library_name, library_strategy, library_selection,
                                 bio_sample, sample_title, taxon_id, taxon_name, age, gender, disease, body_site,
                                 source_name, histological_type, tissue_type, is_tumor, subject_status, hbv_infection,
                                 cirrhosis, tumor_stage, tumor_grade, cell_type, primary_tumor_cirrhosis, sample_type,
                                 clinical_characteristics, phenotype, tumor_type, tumor_status,
                                 cancer_type, diagnosis, health_state, label, description, rna_source
                                 ]

                # disease state missing
                # cause of resections missing

    normal = re.compile(r"([N|n]on-?\_?[T|t]umor)|([N|n]on-?\_?[C|c]ancer)|(normal)|([H|h]ealthy)|([C|c]ontr?ol)"
                        r"|(\_ant)|(adjacent)|([\_|-]N\d)|(\dN;?)|([N|n]or)|([C|c]trl)|((B|(MC))-?\d+[\_]RNA-Seq)",
                        re.IGNORECASE)

    tumor = re.compile(r"([T|t]o?umor)|([H|h]cc)|([H|h]epatocellular)|([C|c]ancer)([D|d]isease)|(carcinoma)|"
                       r"(\_T\d)|(\dT;?)|(icc)|(hgdn)|(lgnd)|(p?dc)|(metastas)", re.IGNORECASE)

    # filter for normal keywords
    for k, v in metadata.items():
        normal_list = [bool(re.search(normal, e)) for e in v]

        if any(normal_list):
            metadata[k].append('normal')
        else:
            metadata[k].append('')

    # filter remaining empty values for tumor keywords
    for k, v in metadata.items():
        tumor_list = [bool(re.search(tumor, e)) for e in v]
        if v[-1] == '' and any(tumor_list):
            v[-1] = 'tumor'

    # create dataframe from rich metadata
    df = pd.DataFrame.from_dict(metadata, orient='index',
                                columns=['Experiment_title', 'Study_ref', 'Library_name', 'Library_strategy',
                                         'Library_selection', 'BioSample', 'Sample_title', 'Taxon_id', 'Taxon_name',
                                         'Age', 'Gender', 'Disease', 'Body_site', 'Source_name', 'Histological_Type',
                                         'Tissue_type', 'Is_tumor', 'Subject_status', 'HBV-infection', 'Cirrhosis',
                                         'Tumor_stage', 'Tumor_grade', 'Cell_type', 'Primary_tumor_cirrhosis',
                                         'Sample_type', 'Clinical_characteristics', 'Phenotype', 'Tumor_type',
                                         'Tumor_status', 'Cancer_type', 'Diagnosis', 'Health_state', 'Label',
                                         'Description', 'Rna_source', 'Condition'
                                         ])

    # save dataframe on disk
    df.reset_index(level=0, inplace=True)
    df.rename({'index': 'RunID'}, axis=1, inplace=True)
    df.to_csv(outpath + '.csv', header=True, index=False)

    # save metadata table additionally as short table
    metadata_df = df[['RunID', 'BioSample', 'Condition', 'Study_ref']].copy()
    metadata_df.rename(
        columns={'RunID': 'FileID', 'BioSample': 'CaseID', 'Condition': 'SampleType', 'Study_ref': 'Project'},
        inplace=True)
    metadata_df.to_csv(outpath + '_short.csv', header=True, index=False)


if __name__ == "__main__":
    sys.exit(main())
