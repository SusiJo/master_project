# This is a script to parse multiple files with TPM values for Differential Gene Expression analysis
# Purpose: enhance performance by using built-in python

import fileUtils.file_handling as fh
import logging
import click
import time
import sys

# Create logger
logger = logging.getLogger('TPM table creator')
# Create console handler
ch = logging.StreamHandler()
# Create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)
logger.setLevel(logging.INFO)


@click.command()
@click.option('-i', '--inpath', prompt='path to stringtie Files',
              help='Path to folder with TPM values',
              required=True)
@click.option('-g', '--genes', prompt='gene ids from rnaseq',
              help='Path to extracted genes from rnaseq file',)
@click.option('-o', '--outpath', prompt='output file',
              help='Output file; Table containing TPM values',
              required=True)
def main(inpath, genes, outpath):
    print("Starting...")
    start_time = time.time()

    # printed to STDOUT on command line
    logger.info('Get stringTie input files')
    allfiles, sample_ids = fh.get_files_and_sample_ids(inpath, '*.txt')

    logger.info('Parse stringTie files')
    gene_tpm, gene_dict = parse_stringTie(allfiles, genes)

    logger.info('Write TPM values to table')
    fh.write_table(gene_dict, sample_ids, gene_tpm, outpath)

    end_time = time.time()
    logger.info('Process finished in ' + str(round(end_time - start_time, 2)) + "sec")


def parse_stringTie(files, genes):
    # Adapted script from Steffen Lemke
    """Read a file and save gene_ids and tmp_values

    :param files: list of files
    :param genes: dict mapping gene_ids to gene_names
    :return: tpm_all_files list of lists
    :return: gene_dict of unique gene_ids to gene_names
    """
    gene_dict = fh.read_unique_genes(genes)

    # initialize dict with keys from gene_dict and empty list as value --> where the tmp values are to be stored
    tpm_all_files = {key: [] for key in gene_dict.keys()}

    for file in files:
        # Read each file
        file = file.split('/')[-1]
        # print(file)
        with open(file, 'r') as f:
            tpm_dict_one_file = {}
            # read file from 2nd line
            for line in f.readlines()[1:]:
                # split lines by tab
                fields = line.strip().split("\t")
                # take "Gene Id" and "Gene Name"
                current_gene_id = fields[0]
                current_tpm = fields[-1]

                # Append gene_id, tmp_value to temporary dict for 1 file
                if current_gene_id not in tpm_dict_one_file.keys():
                    tpm_dict_one_file[current_gene_id] = float(current_tpm)
                else:
                    tpm_dict_one_file[current_gene_id] += float(current_tpm)

            # Append values of 1 file to tmp_values list comprising all files
            for gene_id, tpm_value in tpm_dict_one_file.items():
                tpm_all_files[gene_id].append(tpm_value)

    return tpm_all_files, gene_dict


if __name__ == "__main__":
    sys.exit(main())
