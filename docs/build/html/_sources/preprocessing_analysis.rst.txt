Preprocessing
===============

This chapter outlines how to use the command-line tools.


Parsing stringTie-TPM
*********************

The script merges TPM values from ``rnaseq/stringTieFPKM/`` and maps them to all possible Ensembl gene_ids 
and Ensembl gene_names of the nf-core/rnaseq pipeline.

.. code-block:: bash

    parseTPM.py
    -i/--inpath:  Path to stringtie folder
    -g/--genes:   Path to extracted genes from rnaseq file
    -o/--outpath: Path to output file - tsv table containing merged TPM values

    python parseTPM.py -i <path_to/stringTieFPKM/> -g <path_to/genes_id2name> -o <path_to/mergedTPM>



Parsing FeatureCounts
*********************

Depending on whether the pipeline has run through without issues it is possible to 
directly take the ``rnaseq/featureCounts/merged_gene_counts.txt``.

.. code-block:: bash

    merge_gene_counts.py
    -i/--inpath:      Path to folder with nf-core/rnaseq/featureCounts/merged_gene_counts.txt tables to be merged
    -o/--outpath:     Destination directory. Extension of file filled automatically
    -p/--pseudogenes: Path to file with tab separated file with pseudogenes as EnsemblID

    python merge_gene_counts.py -i <merged_gene_counts.txt> -p <pseudogenes.txt> -o <gene_table.txt>


If the pipeline had errors and stopped, it might be necessary to merge the raw counts
separately from ``rnaseq/featureCounts/gene_counts/``.

.. code-block:: bash

    parseFeatureCounts.py
    -i/--inpath:  	    Path to featureCounts/gene_counts/ 
    -g/--genes:   	    Path to extracted genes from rnaseq file
    --lengths/--no-lengths:  Extract featureLengths (into csv) for raw count to TPM conversion (do this only once) 
    -o/--outpath: 	    Path to output file - tsv table containing merged raw counts 

    python parseFeatureCounts.py -i <path_to/gene_counts/> --no-lengths -g <path_to/genes_id2name> -o <path_to/mergedFeatureCounts>


Merging Replica
***************

If biological replicates are present, it is necessary to merge them into a single observation.
For this purpose the ``merge_replica.py`` script needs the associated metadata table.
It produces two gene_count output tables which can be both input into the clustering analysis.

The outputs are:


* gene_count_table [GeneID, GeneName, ID1, ID2, ...]
* transformed metadata [FileID, CaseID, SampleType, Project]
* combined gene_count_metadata_table [CaseID, SampleType, Condition, Project, UUID, ID1, ID2, ...]


.. code-block:: bash

    merge_replica.py
    -i/--inpath:        Path to merged_gene_counts.txt/ mergedFeatureCounts.txt
    -m/--metadata:      Path to metadata from metadata_processing.py
    -o/--outpath:       Output directory 



Batch Correction
****************

Adjust for known batch effects arising from technical variation with ComBat, ComBat-Seq (svaseq) and ``removeBatchEffect()`` from Limma. 
Needs raw featureCounts as input and converts the counts to TPM units.
The `batch_correction.R <https://github.com/SusiJo/master_project/blob/main/R_notebooks/batch_correction.R>`_ script should be run 
from within the R_docker_container.
Additionally performs Principal Component Analysis for a first overview of the effect of batch correction.

.. code-block:: bash

    # run in R_container

    Rscript container/app/batch_correction.R 
    -c/--counts:    Path to gene_count_table.txt
    -m/--metadata:  Path to metadata.csv
    -l/--lengths:   Path to featureLengths.txt - can be obtained by running parseFeatureCounts.py once with --lengths 
    -t/--title:     Title for plots 
    -o/--outpath:   Destination outpath 



Analysis
=========


Dimension Reduction / Embeddings with PCA, t-SNE and UMAP
*********************************************************

.. note:: 
    The ``mergedTPM/gene_counts.txt`` file has 63,675 unique genes and are highly dimensional. 
    Without pseudogenes there are still 47,224 genes left.
    In order to get a feeling for the structure of the data, it is helpful to 
    reduce the number of dimensions to a lower dimensional embedding. 
    It is advisable to run the script before and after batch correction to get additional plots
    with non-linear dimensionality reduction techniques.


:func:`dim_reduction.py`

* scales the data on the gene features with sklearn.preprocessing.MinMaxScaler() to preserve the variances
* can perform unsupervised linear reduction via Principal Component Analysis - PCA
* can perform unsupervised non-linear reduction via T-distributed Stochastic Neighbor Embedding - t-SNE
* can perform non-linear stochastic reduction with Uniform Manifold Approximation Projection - UMAP
* outputs interactive html plots displaying metadata on hover (produced with plotly) and png images
* can produce a comparative embedding in one image
* can produce a silhouette plot to get an estimate of number of clusters in the data (k-means clustering)

.. code-block:: bash

    dim_reduction.py
    -i/--inpath:      Path to gene_expression_table
    -m/--metadata:    Path to file with metadata for hover_info in plotly plots: [FileID, CaseID, SampleType, Project]
    -o/--outpath:     Path to output images folder
    --pca/--no-pca:   Conditional flag to perform PCA
    --tsne/--no-tsne: Conditional flag to perform t-SNE
    --umap/--no-umap: Conditional flag to perform UMAP
    --comparison/--no-comparison: Conditional flag to plot PCA, t-SNE, UMAP in one plot
    --silhouette/--no-silhouette: Conditional flag to perform silhouette analysis
    -t/--title:                   Title for plots

    python dim_reduction.py -i <mergedTPM> -m <metadata> -o <images/> --pca --tsne --umap --no-comparison --no-silhouette














