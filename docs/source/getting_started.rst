Getting Started
===============

This chapter gives an introduction into the preprocessing steps of RNA-Seq expression data from different public repositories.
It serves as guideline on how download data and metadata, setup a conda environment to use the scripts written in bash, nextflow, python and R.


Features
********

* Download RNA-Seq expression data from repositories
* Convert BAM to FASTQ and use nf-core/rnaseq
* Download metadata from TCGA, ICGC, GTEx, SRA
* Extract metadata into a table in csv format
* Merge TPM values from nf-core/rnaseq/stringTieFPKM
* Merge raw featureCounts from nf-core/rnaseq/featureCounts
* Clustering analyses performing PCA, TSNE, UMAP
* Batch correction with ComBat, CombatSeq



Prerequisites
*************

.. seealso::
   Checkout `qbic-docs <https://pipeline-docs.readthedocs.io/en/latest/index.html>`_

   Assure `nextflow <https://www.nextflow.io/docs/latest/getstarted.html>`_, `docker <https://docs.docker.com/engine/install/centos/>`_, `singularity <https://sylabs.io/guides/3.0/user-guide/installation.html#before-you-begin>`_  are installed



Acquisition of RNA-Seq data
****************************

TCGA
----

Download `GDC transfer tool <https://gdc.cancer.gov/access-data/gdc-data-transfer-tool>`_ to ``/usr/local/bin/`` or another location in ``$PATH``

.. code-block:: bash

   wget 'https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.0_Ubuntu_x64-py3.7_0.zip' -nv -O gdc-client_v1.6.0_Ubuntu_x64-py3.7_0.zip
   unzip gdc-client_v1.6.0_Ubuntu_x64-py3.7_0.zip


Obtain a manifest for data download and an authentication token.
Download data with the following command:

.. code-block:: bash

   gdc-client download -m <manifest> -t <token> -d <outdir>



ICGC
----

.. note:: 
   ICGC currently allows download of data from the Collaboratory Repository withouth having specific cloud access.
   Other repositories need cloud access.

Download `score-client <https://docs.icgc.org/download/guide/#installation-of-the-score-client>`_ (requires OpenJDK11) 

.. code-block:: bash

   wget -O score-client.tar.gz https://artifacts.oicr.on.ca/artifactory/dcc-release/bio/overture/score-client/\[RELEASE\]/score-client-\[RELEASE\]-dist.tar.gz
   tar -xvzf score-client.tar.gz
   cd score-client-5.1.0  # or newer version
   bin/score-client


Get access token and store under ``score-client/conf/application.properties``

Download Data

.. code-block:: bash
  
   score-client --profile collab download --manifest <manifest.tsv> --output-dir <dir>



GTEX / NCBI SRA
---------------

Download data from GTEX project via the Short Read Archive using the `qbic-pipelines/sradownloader <https://github.com/qbic-pipelines/sradownloader>`_.
The current version does not support automatic metadata download yet, but it will be a feature in the future.



Obtain RNA-Seq expression data
******************************

* Run `qbic-pipelines/bamtofastq <https://github.com/qbic-pipelines/bamtofastq>`_  ``v.1.1.0``
* Run `nf-core/rnaseq <https://nf-co.re/rnaseq/1.4.2/usage>`_  ``v.1.4.2``



Metadata 
********

Download `tcga metadata <https://github.com/SusiJo/master_project/blob/main/scripts/tcga_metadata.sh>`_ from **TCGA** via API from `Case endpoint <https://docs.gdc.cancer.gov/API/Users_Guide/Appendix_A_Available_Fields/#case-fields>`_ ``json``

.. code-block:: bash

   bash tcga_metadata.sh <manifest> <json_files>


Download `icgc metadata <https://github.com/SusiJo/master_project/blob/main/scripts/icgc_metadata.sh>`_ from **ICGC** via API from `file and donor endpoint <https://docs.icgc.org/portal/api-endpoints/#!/donors/find>`_ ``json``

.. code-block:: bash
   
   bash icgc_metadata.sh <manifest> <json_file_endpt> <json_donor_endpt>


Download `sra metadata <https://github.com/SusiJo/master_project/blob/main/scripts/metadata.nf>`_ from **NCBI SRA** ``csv, xml``

.. code-block:: bash
   
   nextflow run metadata.nf --run_acc_list <SRA.txt> --outdir <results>


Extract metadata from **TCGA, ICGC, SRA into one table** ``csv`` 

.. note:: 
   Allows multiple paths as input


.. code-block:: bash
  
   metadata_processing.py --icgc <path_to_json> --sra <path_to_csv> --tcga <path_to_json> -o <all_metadata.csv>


.. list-table:: Structure of table output from metadata_processing.py
   :widths: 25 25 25 25
   :header-rows: 1

   * - Sample_ID
     - Case_ID
     - Condition
     - Project_ID
   * - SRR...
     - SAMN...
     - normal
     - GTEX
   * - 6cde....
     - DO...
     - tumor
     - ICGC_PACA_AU


Extract rich metadata from TCGA or ICGC containing additional information on primary diagnosis, 
tumor subtype, gender, vital status, age, survival time, tumor stage, Icd10

.. code-block:: bash

   tcga|icgc_metadata_processing.py -i <inpath> -o <outpath>


Extract rich metadata from NCBI SRA by merging xml and csv 
Option: merge all csv with flag --csv_only

.. code-block:: bash

   csv_metadata_processing.py -i <inpath> --csv_only/--no_csv_only -o <outpath>



Setting-up conda environment
****************************

* Requires ```conda``
* Requires ``python version 3.8.8``

.. code-block:: bash
   
   conda env create -f environment.yml


Activate the environment to run the scripts 

.. code-block:: bash
   
   conda activate python_scripts



Parsing TPM values from
***********************

The script merges TPM values from ``rnaseq/stringTieFPKM/`` and maps them to Ensembl gene_ids and Ensembl gene_names.

.. code-block:: bash

   parseTPM.py
   -i/--inpath:  Path to stringtie folder
   -g/--genes:   Path to extracted genes from rnaseq file
   -o/--outpath: Path to output file - tsv table containing merged TPM values

   python parseTPM.py -i <path_to/stringTieFPKM/> -g <path_to/genes_id2name> -o <path_to/mergedTPM>



Parsing FeatureCounts from
**************************

Depending on whether the pipeline has run through without issues it is possible to 
directly take the ``rnaseq/featureCounts/merged_gene_counts.txt``.

If the pipeline had errors and stopped, it might be necessary to merge the raw counts
separately from ``rnaseq/featureCounts/gene_counts/``.

.. code-block:: bash

   parseFeatureCounts.py
   -i/--inpath:  	    Path to featureCounts/gene_counts/ 
   -g/--genes:   	    Path to extracted genes from rnaseq file
   --lengths/--no-lengths:  Extract featureLengths (into csv) for raw count to TPM conversion (do this only once) 
   -o/--outpath: 	    Path to output file - tsv table containing merged raw counts 

   python parseFeatureCounts.py -i <path_to/gene_counts/> --no-lengths -g <path_to/genes_id2name> -o <path_to/mergedFeatureCounts>



Dimension Reduction/ Clustering analysis
****************************************

.. note:: 
   The ``mergedTPM.txt`` has 63,675 unique genes and is highly dimensional. 
   In order to get a feeling for the structure of the data, it is helpful to 
   reduce the number of dimensions to a lower dimensional embedding. 


:func:`clustering.py`
     

* scales the data on the gene features with sklearn.preprocessing.MinMaxScaler() to preserve the variances
* can perform unsupervised linear reduction via Principal Component Analysis - PCA
* can perform unsupervised non-linear reduction via T-distributed Stochastic Neighbor Embedding - TSNE
* can perform non-linear stochastic reduction with Uniform Manifold Approximation Projection - UMAP
* outputs interactive html plots displaying metadata on hover (produced with plotly) and png images
* can plot a comparative clustering plot in one png image
* can perform silhouette plot analysis to get a good estimate of number of clusters (k-means clustering)

.. code-block:: bash

   clustering.py
   -i/--inpath:      Path to mergedTPM table
   -m/--metadata:    Path to file with metadata for hover_info in plotly plots: [ID, Case_ID, Condition, Project]
   -o/--outpath:     Path to images folder
   --pca/--no-pca:   Conditional flag to perform PCA
   --tsne/--no-tsne: Conditional flag to perform TSNE
   --umap/--no-umap: Conditional flag to perform UMAP
   --comparison/--no-comparison: Conditional flag to plot PCA, TSNE, UMAP in one plot
   --silhouette/--no-silhouette: Conditional flag to perform silhouette analysis

   python clustering.py -i <mergedTPM> -m <metadata> -o <images/> --pca --tsne --umap --no-comparison --no-silhouette



Batch Correction
****************

Adjust for known batch effects arising from technical variation. 

.. warning::
   Currently only tested and implemented in R using svaseq package ``ComBat`` and ``ComBatSeq``
 
`Batch_correction.Rmd <https://github.com/SusiJo/master_project/blob/main/R_notebooks/Batch_correction.Rmd>`_










