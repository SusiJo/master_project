



# master_project

[toc]

## Utilities

* fileUtils contains functions to handle stringTie|featureCounts Files

* `extract_gene_ids2name.py` takes a gtf file as input and extracts all unique geneIDs mapped to their geneNames (i.e. all genes from `nf-core/rnaseq`)

  

**Downloading analysis scripts**

```bash
<score-client|gdc-client> download ... 2>&1 > output.log
```

```bash
./downloading_time.py -i output.log -d <tcga|icgc>
```



## Metadata Scripts

1. **Download json_files TCGA** via API

   ```bash
   ./tcga_metadata.sh <manifest> <json_files_outdir>
   ```

   

2. **Download json_files ICGC**: file and donor endpoint from API

   ```bash
   ./icgc_metadata.sh <manifest> <json_file_endpt> <json_donor_endpt>
   ```

   

3. **Download metadata SRA**: Experiment_info in xml, RunInfo in csv

   ```bash
   nextflow run metadata.nf --run_acc_list <SRA.txt> --outdir <results>
   ```

   

4. **Extract metadata for clustering script**: allows multiple paths as input

   ```bash
   metadata_processing.py -i|--icgc <path_to_json> -s|--sra <path_to_csv> -t|--tcga <path_to_json> -o <all_metadata.csv>
   ```

    

   | Sample_ID | Case_ID | Condition | Project_ID   |
   | --------- | ------- | --------- | ------------ |
   | SRR...    | SAMN... | normal    | GTEX         |
   | 6cde....  | DO...   | tumor     | ICGC_PACA_AU |

   

5.  **Extract rich metadata TCGA|ICGC** 

   ```
   <tcga|icgc>_metadata_processing.py -i <inpath> -o <outpath.csv>
   ```

   **TCGA**

   File_ID, Submitter_ID, Case_ID, Project, Primary_site, Sample_type, Primary_diagnosis, Gender, Vital_status, Age_at_index, Survival_time, Tumor_stage, Icd10

   **ICGC**

   File_ID, File_Name, Donor_ID, Project, Sample_ID, Specimen_ID, Sample_type, Primary_site, Primary_diagnosis, Tumor_subtype, Gender, Vital_status, Age_at_index, Survival_time, Tumor_stage, Icd10



**TODO**: unify 

File_ID, File_Name, Case_ID, Project, Primary_site, Sample_type, Tumor_type, Condition, Gender, Vital_status, Age, Survival_time, Tumor_stage, Icd10



1. **Extract rich metadata SRA** or merge all csv RunInfo files

```bash
csv_metadata_processing.py -i <inpath> --csv_only/--no_csv_only -o <outpath.csv>
```

RunID, Project, BioSample, Library_Strategy, Organism, Gender, Age, Source_name, Subject_status, Tissue, Phenotype, Is_Tumor, StudyDisease, HistologicalType, Condition



## TPM| Feature Counts Scripts

**Set-up venv with python on server**

```bash
python3 -m venv dirname-env
source dirname-env/bin/activate
# install required packages in activated environment
python3 -m pip install -r requirements.txt
```

**Set-up env with Conda on server**

```bash
# from environment file
conda env create -f environment.yml 
# or with specific version
conda create -n snowflakes python=3.8
conda snowflakes activate
```

**Parse TPM**

```bash
python parseTPM.py -i stringTieFPKM/ -g genes_id2name.txt -o mergedTPM.txt
```

**Parse FeatureCounts**

```bash
# :param lengths|no-lengths: boolean flag: if set, extracts gene-lengths from feature counts
python parseFeatureCounts.py -i featureCounts/ -g genes_id2name.txt -o mergedFeatureCounts.txt
```

**MergeReplica**

```bash

```



## Preprocessing Scripts

**Clustering script**

* performs PCA, TSNE, UMAP, ComparativePlot, SilhouettePlot

```bash
python clustering.py -i <mergedTPM.txt> -m <metadata.csv> -o <outpath_images> 
[ -pca, -tsne, -umap, -comparison, -silhouette ]
```



## Batch removal

**R packages**

* Combat, CombatSeq
* removeBatchEffect

```R

```

**Python packages**

* Combat.py
* Pycombat.py

 



## ML jupyter notebooks

**Classification tumor vs. Normal**

