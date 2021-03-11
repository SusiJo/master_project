#!/bin/bash

# Purpose: download metadata from TCGA using the API 
# $1 = path to manifest
# $2 = path to output folder

runinfo='https://api.gdc.cancer.gov/files/UUID?expand=cases,cases.demographic,cases.samples,annotations,cases.diagnoses,cases.diagnoses.treatments,cases.project,cases.project.program,analysis.metadata.read_groups&pretty=true'

sed 1d $1 | cut -f1 > tcga_ids.txt

while IFS=, read -r line
do 
    url=$(echo ""$runinfo"" | sed "s/UUID/$line/")
    curl $url -s --output $2/$line.json
done < tcga_ids.txt

echo Done!
rm tcga_ids.txt
