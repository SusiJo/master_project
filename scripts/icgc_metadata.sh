#!/bin/bash

# Purpose: download metadata from ICGC using the API 
# $1 = path to manifest
# $2 = path to output folder with file_endpoint
# $3 = path to output folder with donor_endpoint

donorinfo='https://dcc.icgc.org/api/v1/donors/DONOR?include=specimen'
metadata='https://dcc.icgc.org/api/v1/repository/files/FILE'

# skip first line with header
sed 1d $1 | cut -f2,5,9 | sed '/TopHat/d' > icgc_ids.txt     # in some manifests donor is at field 8 (not 9) 
while IFS=, read -r line
do 
     # Manifest Column2 = FileID
     fileid=$(echo "$line"| awk -F" "  '{print $1}')
     
     # Manifest Column5 = FileName.ext
     # Substring depends on naming convention of project
     # PACA-CA example (first 32 letters)
     # 4a66f8679e9b414593273824c25e7d5c.SWID_280567_MPCC_0006_Pa_C_PE_250_MR_120925_SN1068_0099_AD0V56ACXX_NoIndex_L004_R2_001.fastq.gz
     # donortoreplace=$(echo "$line"| awk -F" "  '{print $3}')
     # filename=$(echo "$line"| awk -F" "  '{print substr($2, 1, 32)}')
     # PACA-AU example (7-36)
     # PCAWG.e03743ed-526c-4f13-81f4-0f32faca2110.STAR.v1.bam
    
 
     # Get File Name and Donor
     filename=$(echo "$line"| awk -F" "  '{print substr($2, 7, 36)}')
     donortoreplace=$(echo "$line"| awk -F" "  '{print $3}')
     
     # Download metadata FILE ENDPOINT  
     mkdir -p $2
     url=$(echo ""$metadata"" | sed "s|FILE|$fileid|g")
     curl $url -s --output $2/$fileid.$filename.json
	 
     # Download metadata DONOR ENDPOINT
     mkdir -p $3
     donorurl=$(echo ""$donorinfo"" | sed "s|DONOR|$donortoreplace|g")
     curl $donorurl -s --output $3/$fileid.$filename.$donortoreplace.json
done < icgc_ids.txt

echo Done!
rm icgc_ids.txt

# Manifest ICGC
# repo_code, file_id, object_id, file_format, file_name, file_size, md5sum, index_object_id, donor_id/donor_count, 
# project_id/project_count, study
