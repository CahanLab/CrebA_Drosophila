#!/bin/bash
output_path=../output/bound_genes_regulatory_regions/up_genes/fly_factor_survey
mkdir -p ${output_path}

cat ../input/motif_databases.12.23/FLY/fly_factor_survey.meme | grep MOTIF > ${output_path}/motif_database.txt 
while read p; do 
    motif_id=$(echo ${p} | cut -d ' ' -f2)
    motif_name=$(echo ${p} | cut -d ' ' -f3)
    fimo --o ${output_path}/${motif_id}---${motif_name} --motif ${motif_id} ../input/motif_databases.12.23/FLY/fly_factor_survey.meme ../output/bound_genes_regulatory_regions/up_genes/bound_regulatory_regions_seq.fasta
done < ${output_path}/motif_database.txt