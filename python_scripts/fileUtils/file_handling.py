
import os.path
import glob


def get_files(inpath, ext):
    """Create list of files

    :param inpath: to folder
    :param ext: of files
    :return: files: list
    """
    os.chdir(inpath)
    files = [os.path.abspath(os.path.basename(f)) for f in glob.glob(inpath + ext)]
    return files


def get_files_and_sample_ids(inpath, ext):
    """Create list of files

    :param inpath: to files
    :param ext: file extension
    :return: files: as list
    :return: sample_ids: as list
    """
    os.chdir(inpath)
    files = [os.path.abspath(os.path.basename(f)) for f in glob.glob(inpath + ext)]
    sample_ids = []
    for i, file in enumerate(files):
        file = file.split('/')[-1]
        if file.startswith('PCAWG'):
            # print(file.split(".")[1])
            sample_ids.append(file.split(".")[1])
        else:
            # print(file.split("_")[0])
            sample_ids.append(file.split("_")[0])

    return files, sample_ids


def read_unique_genes(genes):
    """Read genes from rnaseq pipeline

    :param genes: genes_id mapped to name
    :return: gene_dict: mapping ids to names
    """
    gene_dict = {}
    with open(genes, 'r') as f:
        for line in f.readlines():
            fields = line.strip().split("\t")
            gene_id = fields[0]
            gene_name = fields[1]
            gene_dict[gene_id] = gene_name
    return gene_dict


def write_table(gene_dict, sample_ids, all_files, outpath):
    """Create tab separated table with merged values, adapted from Steffen Lemke

    :param gene_dict: gene_ids mapped to gene_names
    :param sample_ids: array
    :param all_files: dict
    :param outpath: path where to store outfile, tsv
    """
    with open(outpath, 'w') as table:
        # .join(list) iterates over ids
        table.write("GeneID\tGeneName\t" + "\t".join(sample_ids) + "\n")

        all_gene_ids = list(all_files.keys())
        all_gene_ids.sort()

        for gene_id in all_gene_ids:
            gene_name = gene_dict[gene_id]
            table.write(gene_id + "\t" + gene_name + "\t")
            tpm_values_gene_id_str = [str(x) for x in all_files[gene_id]]
            table.write("\t".join(tpm_values_gene_id_str) + "\n")


def read_tpm(infile):
    """Reads TPM values from merged file into floats

    :param infile: merged TPM values
    :return: sample_ids
    :return: gene_ids
    :return: gene_names
    :return: data_list
    """
    tmp = []
    gene_names = []
    gene_ids = []
    header = []
    with open(infile, 'r') as file:

        for i, line in enumerate(file):
            fields = line.strip().split("\t")
            if i == 0:
                # first line: skip first two fields "Gene ID", "Gene Name"
                header.append(fields[2:])
            if i > 0:
                # other lines: create lists
                current_gene = fields[0]
                current_gene_name = fields[1]
                tpm_values = fields[2:]
                # in case that a gene has no tmp values --> skip that line! (alternative: fill with zeros)
                if len(tpm_values) == 0:
                    continue
                gene_ids.append(current_gene)
                gene_names.append(current_gene_name)
                tmp.append(tpm_values)

        # header contains sample ids
        sample_ids = [val for sublist in header for val in sublist]

    file.close()

    # convert nested list elements to floats
    data_list = [list(map(float, item)) for item in tmp]

    return sample_ids, gene_ids, gene_names, data_list


def parse_csv(files):
    """Parse csv files

    :param files: list of files
    :return: srr_metadata: metadata dictionary
    """
    srr_metadata = {}
    is_tumor = ''
    for file in files:
        # Read each file
        file = file.split('/')[-1]
        with open(file, 'r') as f:
            # strip lines to not read empty  lines
            lines = list(line for line in (l.strip() for l in f) if line)
            for i, l in enumerate(lines):
                # skip header
                if i == 0:
                    continue
                fields = l.split(",")
                run = fields[0]
                bio_project = fields[21]
                bio_sample = fields[25]
                if fields[36] == 'no':
                    is_tumor = 'normal'
                if fields[36] == 'yes':
                    is_tumor = 'tumor'
                if run not in srr_metadata:
                    srr_metadata[run] = [bio_sample, is_tumor, bio_project]

    return srr_metadata
