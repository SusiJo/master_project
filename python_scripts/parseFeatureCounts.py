# This is a script to parse multiple files from featureCounts to perform batch correction

import click
import logging
import time
import sys
import pandas as pd
from pathlib import Path
import fileUtils.file_handling as fh

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
@click.option('-i', '--inpath', prompt='path to count Files',
              help='Path to folder with featureCount values',
              required=True)
@click.option('-g', '--genes', prompt='gene ids from rnaseq',
              help='Path to extracted genes from rnaseq file', )
@click.option('-o', '--outpath', prompt='output file',
              help='Output file; Table containing featureCount values',
              required=True)
@click.option('--lengths/--no-lengths', prompt='extract lengths of genes',
              help='writes table with gene lengths if true, default=False',
              default=False)
def main(inpath, genes, outpath, lengths):
    print("Starting...")
    start_time = time.time()

    # printed to STDOUT on command line
    logger.info('Get featureCount input files')
    allfiles, sample_ids = fh.get_files_and_sample_ids(inpath, '*.txt')

    logger.info('Parse featureCount files')
    gene_counts, gene_dict, gene_lengths = parse_featureCounts(allfiles, genes)

    logger.info('Write count values to table')
    table = fh.write_table(gene_dict, sample_ids, gene_counts, outpath)

    if lengths:
        write_lengths(gene_lengths, outpath)

    end_time = time.time()
    logger.info('Process finished in ' + str(round(end_time - start_time, 2)) + "sec")


def parse_featureCounts(files, genes):
    """Read a file and the sample ID and save gene_ids and tmp_values, adapted from Steffen Lemke

    :param files: list of files
    :param genes: dict mapping gene_ids to gene_names
    :return: counts_all_files
    :return: gene_dict
    :return: geneLengths_all_files : dict mapping ids to lengths
    """
    gene_dict = fh.read_unique_genes(genes)

    # initialize dict with keys from gene_dict and empty list as value --> where the tmp values are to be stored
    counts_all_files = {key: [] for key in gene_dict.keys()}
    geneLengths_all_files = {key: [] for key in gene_dict.keys()}

    for file in files:
        file = file.split('/')[-1]
        # Read each file
        with open(file, 'r') as f:
            count_dict_one_file = {}
            # read file from 2nd line
            for line in f.readlines()[2:]:
                # split lines by tab
                fields = line.strip().split("\t")
                # take "Gene Id" and "Gene Name"
                current_gene_id = fields[0]
                current_length = fields[5]
                # read count at last column
                current_count = fields[-1]

                # Append gene_id, tmp_value to temporary dict for 1 file
                if current_gene_id not in count_dict_one_file.keys():
                    count_dict_one_file[current_gene_id] = int(current_count)
                    # save exon gene lengths: they are identical for all featureCounts outputs
                    geneLengths_all_files[current_gene_id] = current_length
                else:
                    count_dict_one_file[current_gene_id] += int(current_count)

            # Append values of 1 file to tmp_values list comprising all files
            for gene_id, tpm_value in count_dict_one_file.items():
                counts_all_files[gene_id].append(tpm_value)

    return counts_all_files, gene_dict, geneLengths_all_files


def write_lengths(geneLengths_all_files, outpath):
    """Create mini-table with gene-lengths from RNAseq pipeline from FeatureCounts output
        Comprises all non-overlapping bases in exons belonging to the same gene

    :param geneLengths_all_files:  gene_ids mapped to gene_lengths
    :param outpath: path where to store outfile, txt
    """
    table = pd.DataFrame.from_dict(geneLengths_all_files, orient='index', columns=['GeneLength'])
    table.reset_index(level=0, inplace=True)
    table.rename({'index': 'GeneID'}, axis=1, inplace=True)
    p = Path(outpath)
    ext = p.suffix
    new_outpath = p.rename(Path(p.parent, "featureCounts_geneLengths" + ext))
    # print(new_outpath)
    table.to_csv(new_outpath, index=False)


if __name__ == "__main__":
    sys.exit(main())
