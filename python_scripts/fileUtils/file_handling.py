
import os.path
import glob


def get_files(inpath, ext):
    """Create list of files
    :param inpath
    :param ext of files
    @:return: list of files"""
    os.chdir(inpath)
    files = [os.path.abspath(os.path.basename(f)) for f in glob.glob(inpath + ext)]
    return files


def get_files_and_sample_ids(inpath, ext):
    """Create list of files
    @:param: inpath to files
    @:param: file extension
    @:return list of files, list of sample_ids"""
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
    # print(sample_ids)

    return files, sample_ids


def read_unique_genes(genes):
    """Read genes from rnaseq pipeline
    @:param: genes
    @:return: gene_dict: mapping ids to names"""
    gene_dict = {}
    with open(genes, 'r') as f:
        for line in f.readlines():
            fields = line.strip().split("\t")
            gene_id = fields[0]
            gene_name = fields[1]
            gene_dict[gene_id] = gene_name
    return gene_dict


def write_table(gene_dict, sample_ids, all_files, outpath):
    """Create tab separated table with merged values
    From Steffen Lemke
    @:param: gene_dict: gene_ids mapped to gene_names
    @:param: sample_ids
    @:param: allfiles
    @:param: outpath"""

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
    @:param: infile: merged TPM values
    @:return: sample_ids, gene_ids, gene_names, data_list
    """

    tmp = []
    gene_names = []
    gene_ids = []
    header = []
    sample_ids = []
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


