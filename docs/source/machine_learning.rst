Machine-Learning on RNA-Seq data
================================

This chapter describes how to use the machine learning script and use different KEGG genes as filter.


Machine Learning
****************

This script should be run within the docker container using the ``python_scrips/environment.yml``.
For efficient investigation of different algorithms a wrapper bash command is outlined here.
It executes the machine_learning script for all implemented algorithms and writes STDOUT and STDERR to a log file.

Outputs: 

* Plot of 10 most important features
* Feature importances as csv file
* Plot of confusion matrix
* Plot of ROC-AUC curve

.. code-block:: bash

    # Bash wrapper script

    declare -a classifiers=("LinearSVC" "SVC" "RandomForest" "MultiLayerPerceptron") 

    for val in ${classifiers[@]}; do

        echo Starting to train $val classifier

        python container/app/machine_learning_tool.py \
        -i/--inpath <gene_counts/TPM.txt> \
        -m /--metadata <metadata.csv> \
        -a/--algorithm $val \
        [-k/--kegg <KEGG_filter.txt>] \
        -o/--outpath <destination directory> \
        -c/--cores <INT>  \
        -t/-title <tissue_dataset> 2>&1 >> out.log 

    done 



KEGG database
*************

For obtaining human or cancer genes from the KEGG database the following sources can be used.
The tables have to be merged in a database-like fashion to be input in the ML script.

* Source of human KEGG [genes](http://rest.kegg.jp/list/hsa)
* Source of human KEGG [pathways](http://rest.kegg.jp/list/pathway/hsa)
* Source of NCBI Entrez Identifiers to KEGG [identifiers](http://rest.kegg.jp/conv/hsa/ncbi-geneid)
* Source of mapping NCBI Entrez to Ensembl [ids](https://www.genenames.org/cgi-bin/download/custom?col=gd_pub_eg_id&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit)

- 05200 [Pathways in cancer](https://www.kegg.jp/kegg-bin/show_pathway?hsa05200)
- 05202 [Transcriptional misregulation in cancer](https://www.kegg.jp/kegg-bin/show_pathway?hsa05202)
- 05206 [MicroRNAs in cancer](https://www.kegg.jp/kegg-bin/show_pathway?hsa05206) 
- 05205 [Proteoglycans in cancer](https://www.kegg.jp/kegg-bin/show_pathway?hsa05205)
- 05204 [Chemical carcinogenesis](https://www.kegg.jp/kegg-bin/show_pathway?hsa05204)
- 05203 [Viral carcinogenesis](https://www.kegg.jp/kegg-bin/show_pathway?hsa05203)
- 05230 [Central carbon metabolism in cancer](https://www.kegg.jp/kegg-bin/show_pathway?hsa05230)
- 05231 [Choline metabolism in cancer](https://www.kegg.jp/kegg-bin/show_pathway?hsa05231)
- 05235 [PD-L1 expression and PD-1 checkpoint pathway in cancer](https://www.kegg.jp/kegg-bin/show_pathway?hsa05235)



.. list-table:: Excerpt of a KEGG_filter.txt file
    :widths: 20 20 60
    :header-rows: 1 

    * - ensembl_id
      - GeneName
      - pathway
    * - ENSG00000000419	
      - DPM1	
      - ['path:hsa01100', 'path:hsa00510']
    * - ENSG00000000938	
      - FGR	
      - ['path:hsa04062']
    * - ENSG00000000971	
      - CFH	
      - ['path:hsa04610', 'path:hsa05150']
    * - ENSG00000001036	
      - FUCA2	
      - ['path:hsa04142', 'path:hsa00511']



Feature importances
*******************

Annotation of cancer prognostics can be done with the ``get_protein_atlas.py`` script.
It deploys programmatic data access to the API of `The Human Protein Atlas <https://www.proteinatlas.org/>_` and retrieves information in JSON format.

.. code-block:: bash

    get_protein_atlas.py
    -i/--inpath:            Inpath to table with feature importance values (CSV)
    -c/--cancer:            Cancer type used for query
    -o/--outpath:           Outpath for annotated feature importance table (CSV) 
    --latex:                If set, print annotated table in latex format to STDOUT

    python get_protein_atlas.py -i <feature_importance.csv> -c <Pancreatic/Liver>  -o <annotated_feature_importance.csv>