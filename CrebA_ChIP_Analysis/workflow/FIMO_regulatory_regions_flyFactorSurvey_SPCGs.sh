#!/bin/bash
mkdir -p ../output/SPCG_regulatory_regions/fly_factor_survey/
cat ../input/motif_databases.12.23/FLY/fly_factor_survey.meme | grep MOTIF > ../output/SPCG_regulatory_regions/fly_factor_survey/motif_database.txt 
while read p; do 
    motif_id=$(echo ${p} | cut -d ' ' -f2)
    motif_name=$(echo ${p} | cut -d ' ' -f3)
    fimo --o ../output/SPCG_regulatory_regions/fly_factor_survey/${motif_id}---${motif_name} --motif ${motif_id} ../input/motif_databases.12.23/FLY/fly_factor_survey.meme ../output/SPCG_regulatory_regions/spcg_regulatory_regions_seq.fasta
done < ../output/SPCG_regulatory_regions/fly_factor_survey/motif_database.txt