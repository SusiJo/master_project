"""Script to programmatically access The Human Protein Atlas for Ensembl Gene IDs
    and to obtain pathology prognostics for cancer
"""

import requests
import click
import logging
from time import time
import json
import pandas as pd

# Create logger
logger = logging.getLogger('Protein Atlas')
# Create console handler
ch = logging.StreamHandler()
# Create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add ch to logger
logger.addHandler(ch)
logger.setLevel(logging.INFO)


@click.command()
@click.option('-i', '--inpath', prompt='path to feature_importance table (CSV)',
                help='Path to feature_importance table',
                required=True)
@click.option('-c', '--cancer', help='Enter cancer type to be searched for without quotes',
                type=click.Choice(['Pancreatic', 'Liver'], case_sensitive=True))
@click.option('-o', '--outpath', prompt='path to output table',
                help='Path for annotated feature_importance table. Extension filled automatically (CSV)',
                required=True)
@click.option('--latex', help='If set, prints tables to STDOUT in latex format', is_flag=True)
def main(inpath, cancer, outpath, latex):
    start_time = time()

    logger.info('Parse feature_importance table...')
    # read feature-importance table
    try:
        df = pd.read_csv(inpath)
        # read columns into lists
        feature_importance = df.iloc[:50, 0].tolist()
        gene_ids = df.iloc[:50, 1].tolist()
        gene_name = df.iloc[:50, 2].tolist()

        logger.info('Retrieving attributes...')
        attributes = [get_attribute(i, cancer) for i in gene_ids]
        prognostic_type = [i[0] for i in attributes]
        prognostics = [i[1] for i in attributes]

        # create new dataframe
        df_top50 = pd.DataFrame({'Feature_Importance': feature_importance, 'GeneID': gene_ids,
                                    'GeneName': gene_name, 'Prognostic_Type': prognostic_type,
                                    'Is_prognostic': prognostics})

        logger.info('Write output table...')
        df_top50.to_csv(outpath + "top_50.csv", index=False)

        # select only where prognostic is true
        df_final = df_top50[df_top50['Is_prognostic'] == True]
        df_final.to_csv(outpath + "_prognostic.csv", index=False)

        # print tables in latex format
        if latex:
            print(df_top50.to_latex(index=False))
            print(df_final.to_latex(index=False))

        end_time = time()
        logger.info('Process finished in ' + str(round(end_time - start_time, 2)) + "sec")

    except OSError:
        print('Cannot find file on disk. Please try another path.')


def get_attribute(ensembl_id, cancer):
    """Get prognostic attributes from The Human Protein Atlas via API

    :param ensembl_id: list
    :param cancer: one of Pancreas, Liver
    :return: tuple of [prognostic_type (str), is_prognostic (boolean)]
    """
    endpt = "https://www.proteinatlas.org/{}.json".format(ensembl_id)
    response = requests.post(endpt, headers={"Content-Type": "application/json"})
    try:
        data = json.loads(response.text)
        cancer_type = 'Pathology prognostics - ' + cancer + ' cancer'

        prognostic_type = data[cancer_type]['prognostic type']
        prognostic = data[cancer_type]['is_prognostic']

    except Exception:
        prognostic_type = 'NA'
        prognostic = 'NA'

    return [prognostic_type, prognostic]


if __name__ == "__main__":
    main()
