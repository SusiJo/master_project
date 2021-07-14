# Purpose of this script is to merge merged_gene_counts tables from featureCounts

import click
import logging
import time
import sys
import pandas as pd
import re
import fileUtils.file_handling as fh
from functools import reduce

# Create logger
logger = logging.getLogger('FeatureCounts table creator')
# Create console handler
ch = logging.StreamHandler()
# Create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)
logger.setLevel(logging.INFO)


@click.command()
@click.option('-i', '--inpath', prompt='Please input path to folder with files to be merged',
              help='Please specify input path to folder with files to be merged',
              required=True)
@click.option('-o', '--outpath', prompt='Please enter an output path without extension',
              help='Please enter an output path for merged featureCounts table. Extension is filled automatically.',
              required=True)
@click.option('-p', '--pseudogenes',
              help='Please provide list of pseudogenes if they should be filtered out', required=False)
def main(inpath, outpath, pseudogenes):

    start_time = time.time()

    logger.info('Getting files')
    allfiles = fh.get_files(inpath, '*.txt')

    logger.info('Merging files')
    if pseudogenes:
        pseudogenes_df = pd.read_csv(pseudogenes, sep='\t', header=None)
        pseudogene_list = pseudogenes_df.iloc[:, 0]
        print("Number of pseudogenes: ", len(pseudogene_list))
        merge_frames(allfiles, outpath, pseudogene_list=pseudogene_list)
    else:
        merge_frames(allfiles, outpath, pseudogene_list=None)

    end_time = time.time()
    logger.info('Process finished in ' + str(round(end_time - start_time, 2)) + "sec")


def merge_frames(files, outpath, pseudogene_list):
    """
    Function to glue merged_gene_counts.txt from featureCounts together
    Tables should be in the same folder

    :param files: list of files
    :param outpath: enter filename for output
    :param pseudogene_list: enter a list of pseudogenes to be filtered out
    :return: DataFrame with merged featureCounts values, tab-separated, txt
    """

    dfs = list()
    for filename in files:
        # print(filename)
        df = pd.read_csv(filename, sep='\t')
        df.rename(columns={"Geneid": "GeneID", "gene_name": "GeneName"}, inplace=True)
        dfs.append(df)

    # merge dataframes column-wise based on 2 key columns, axis=1
    df_final = reduce(lambda left, right: pd.merge(left, right, on=['GeneID', 'GeneName']), dfs)
    df_final.sort_values(by='GeneID', inplace=True)

    # if duplicates present, remove them
    long_sample_ids = df_final.columns[2:].tolist()
    mask = [i for i in long_sample_ids if i.endswith('_y')]
    to_select = df_final.columns[~df_final.columns.isin(mask)]
    # remove duplicates
    df_no_dups = df_final[to_select].copy()
    sample_ids_no_dups = df_no_dups.columns[2:].tolist()

    # make ids shorter
    col_names = ['GeneID', 'GeneName']
    short_id_list = [i.split(".")[1] if i.startswith('PCAWG') else re.split("((_1)?Aligned)|(.sra)|(_gdc)", i)[0]
                     for i in sample_ids_no_dups]

    col_names.extend(short_id_list)
    df_no_dups.columns = col_names

    # save dataframe as output
    outpath_suffix = outpath + '.txt'
    df_no_dups.to_csv(outpath_suffix, sep='\t', index=False)

    # apply filter to remove pseudogenes if list provided and save additionally
    if pseudogene_list is not None:
        inverse_boolean_series = ~df_no_dups.GeneID.isin(pseudogene_list)
        filtered_df = df_no_dups[inverse_boolean_series].copy()
        filtered_df.sort_values(by='GeneID', inplace=True)
        filtered_outpath = outpath + "_reduced_geneset.txt"
        filtered_df.to_csv(filtered_outpath, sep='\t', index=False)

    # Print INFO
    print("Number of samples merged_gene_table: ", df_no_dups.columns[2:])


if __name__ == "__main__":
    main()
