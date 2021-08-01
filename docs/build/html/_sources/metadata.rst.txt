Metadata acquisition
====================


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
    ICGC currently allows download of data from the Collaboratory Repository without having specific cloud access.
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

   * - FileID
     - CaseID
     - SampleType
     - Project
   * - FileID1
     - SampleID1
     - normal
     - Project1
   * - FileID1
     - SampleID2
     - tumor
     - Project2


Extract rich metadata from TCGA or ICGC containing additional information on primary diagnosis, 
tumor subtype, gender, vital status, age, survival time, tumor stage, Icd10

.. code-block:: bash

    tcga|icgc_metadata_processing.py -i <inpath> -o <outpath>


Extract rich metadata from NCBI SRA xml file and annotate with conditions.

.. code-block:: bash

    xml_soup.py -x <inpath_xml_dir> -o <outpath>
