""" Script to merge stringTieFPKM TPM or featureCounts raw count values
    outputs: - merged_replica_table
             - merged_gene_table
             - new_metadata with new ids
             - new_metadata with batch information encoded as INT
"""

# imports
import click
import sys
import logging
import time
import pandas as pd

# Create logger
logger = logging.getLogger('Replica merger')
# Create console handler
ch = logging.StreamHandler()
# Create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)
logger.setLevel(logging.INFO)


@click.command()
@click.option('-i', '--inpath', prompt='input path table with replica',
              help='input path to merged TPM/featureCounts table with replica',
              required=True)
@click.option('-m', '--metadata', prompt='path to metadata files',
              help='Path to table in csv format',
              required=True)
@click.option('-o', '--outpath', prompt='output table',
              help='Output table with median merged replica (TXT). Extension filled automatically',
              required=True)
def main(inpath, metadata, outpath):
    start_time = time.time()

    logger.info('Reading metadata...')
    meta_table, meta_dict = read_metadata(metadata)

    logger.info('Writing output to table...')
    merge_replica(inpath, meta_table, outpath)

    end_time = time.time()
    logger.info('Process finished in ' + str(round(end_time - start_time, 2)) + "sec")


def read_metadata(meta):
    """
    Function to read metadata table and connect case_id to file_id
    Format:
    FileID, CaseID, Sample_Type, Project

    :param: meta: inpath to metadata table
    :return: metadata: as dataFrame
    :return: dictionary: mapping case_id to sample_ids
    """

    metadata = pd.read_csv(meta, sep=",")
    file_id = metadata.iloc[:, 0]
    case_id = metadata.iloc[:, 2]
    id_dict = dict(zip(file_id, case_id))

    # dict mapping {BioSample: (SRR1, SRR2), ...}
    new_dict = {}
    for k, v in id_dict.items():
        if v not in new_dict.keys():
            new_dict[v] = [k]
        else:
            new_dict[v].append(k)

    return metadata, new_dict


def merge_replica(inpath, meta_table, outpath):
    """
    Function to merge replica with same BioSample = Case_ID and condition

    :param inpath: path to TPM/featureCounts merged table
    :param meta_table: metadata table, dataFrame
    :param outpath: path where to store merged_file, csv
    :return: merged_replica_file.csv
    """
    pd.set_option("display.max_rows", 20, "display.max_columns", 10)

    # read tpm values
    df = pd.read_csv(inpath, sep="\t")
    gene_name = list(df.iloc[:, 1])
    # print(gene_name)
    # gene_ids = df.index
    df = df.drop(['GeneName'], axis=1)
    # transpose df
    dft = df.T
    sample_ids = dft.index[1:].tolist() # = sorter
    # print("sorter\n ", sample_ids)
    new_header = dft.iloc[0]
    dft = dft[1:]
    dft.columns = new_header

    # print(dft.iloc[:10, :2])

    # print("Meta before sorting\n", meta_table)
    # sort metadata according to sample_ids
    meta_table.FileID = meta_table.FileID.astype("category")
    meta_table.FileID.cat.set_categories(sample_ids, inplace=True)
    meta_table = meta_table.sort_values('FileID')
    # print("Meta after sorting\n", meta_table)

    # attach info to merge replica on
    file_id = list(meta_table.iloc[:, 0])
    case_id = list(meta_table.iloc[:, 1])
    sample_type = list(meta_table.iloc[:, 2])
    project = list(meta_table.iloc[:, 3])
    p_dict = dict(zip(case_id, project))

    mapping = zip(sample_ids, zip(file_id, case_id))

    dft.insert(0, 'CaseID', case_id)
    dft.insert(1, 'SampleType', sample_type)

    dft.reset_index(drop=True, inplace=True)
    dft = dft.rename_axis(None, axis=1)     # remove index name
    dft.insert(0, 'FileID', file_id)

    # make Multi-Index for merging replica based on case_id and sample_type
    dft.set_index(['FileID', 'CaseID', 'SampleType'], inplace=True)

    # convert object type numbers to numeric
    dft = dft.apply(pd.to_numeric)

    # merge replica based sample_id and condition
    df_final = dft.groupby(['CaseID', 'SampleType'], axis=0).median().reset_index()

    # attach condition info as numbers
    condition_dict = {'normal': 0, 'tumor': 1}
    df_final.insert(2, 'Condition', df_final['SampleType'].map(condition_dict))
    # attach project info
    df_final.insert(3, 'Project', df_final['CaseID'].map(p_dict))

    # make new unique ids
    unique_ids = ['id'+str(i) for i in range(len(df_final))]
    df_final.insert(4, 'UniqueID', unique_ids)
    df_final.fillna(0, inplace=True)

    outpath_suffix = outpath + '.txt'
    df_final.to_csv(outpath_suffix, sep="\t", index=False)

    # make dataframe look like before for clustering + batch correction
    df_new = df_final.drop(['CaseID', 'SampleType', 'Condition', 'Project'], axis=1)
    dft_new = df_new.T

    dft_new.drop('UniqueID', inplace=True)
    dft_new.index.set_names('GeneID', inplace=True)
    dft_new.columns = unique_ids
    dft_new.insert(0, 'GeneName', gene_name)

    dft_new.to_csv(outpath+"_gene_table.txt", sep="\t")

    # create new metadata table with new ids
    new_case_ids = list(df_final.CaseID)
    new_sample_types = list(df_final.SampleType)
    new_project = list(df_final.Project)

    metadata_new = pd.DataFrame({'FileID': unique_ids, 'CaseID': new_case_ids,
                                 'SampleType': new_sample_types, 'Project': new_project})

    metadata_new.to_csv(outpath + "_metadata_new.csv", index=False)

    # annotate metadata with Batches
    project_series = pd.Series(metadata_new.Project)
    projects_set = set(metadata_new.Project)
    batches = [i for i in range(1, len(projects_set) + 1)]
    project_batch_dict = dict(zip(projects_set, batches))
    metadata_new['Batch'] = project_series.map(project_batch_dict)
    metadata_new.to_csv(outpath + "_metadata_new_batches.csv", sep=';', index=False)

    # print info
    print("number of samples gene_count_table ", len(dft_new.columns[2:]))
    print("number of samples in metadata table ", len(metadata_new))
    print("No. Normal samples: ", len(metadata_new[metadata_new.SampleType == 'normal']))
    print("No. Tumor samples: ", len(metadata_new[metadata_new.SampleType == 'tumor']))


if __name__ == "__main__":
    main()

