# Purpose of this script is to read all genes that can be found by the nf-core/rnaseq pipeline
# into a list with Ensembl IDs and gene names

# imports

import re


def lines(genes):
    """Generator for lines in gtf file"""
    with open(genes, 'r') as g:
        for line in g:
            yield parse(line)


def parse(line):
    """Parse ;-sep line of gtf file"""
    gene_dict = {}
    # split by tab and only save last column with attribute list
    info_list = line.strip().split("\t")[-1]
    # info format: gene_name "DDX11L1"; gene_id "ENS....";
    # split by ;
    attributes = [x for x in re.split(";", info_list) if x.strip()]

    for i, attribute in enumerate(attributes):
        # split by whitespace --> key = gene_name : value = "DDX11L1"
        key, value = attribute.split()
        # strip trailing characters
        key = key.strip()
        # strip of "
        value = value.strip('"')
        gene_dict[key] = value

    return gene_dict


def create_dict(path):
    """Save unique gene ids and names to dict"""
    id_name_dict = {}
    # iterate through generator object
    gene_dict = lines(path)
    for line in gene_dict:
        # value of 'gene_id' becomes key in new dict
        key = line.get('gene_id')
        value = line.get('gene_name')
        if key not in id_name_dict:
            # dict[key] = value
            # 'ENSG00000132781' : 'MUTYH'
            id_name_dict[key] = value

    return id_name_dict


def write_genes(path, outfile):
    """Call other functions"""
    gene_id_name_dict = create_dict(path)
    """Write dictionary to txt file"""
    with open(outfile, 'w') as file:
        for key, value in gene_id_name_dict.items():
            file.write("{0}{1}{2}{3}".format(key, "\t", value, "\n"))

    # sorted_dict = {}
    # for key in sorted(gene_id_name_dict.keys()):
    #     sorted_dict[key] = gene_id_name_dict[key]


if __name__ == "__main__":
    # print generator in python3
    # print(next(info))

    infile = "/Users/susanne/Documents/Master/genes.gtf"
    write_genes(infile, '/Users/susanne/Documents/Master/genes_id2name.txt')





