#!/bin/bash
mkdir -p ../output/SPCG_promoters/fly_factor_survey/
cat ../input/motif_databases.12.23/FLY/fly_factor_survey.meme | grep MOTIF > ../output/SPCG_promoters/fly_factor_survey/motif_database.txt 
while read p; do 
    motif_id=$(echo ${p} | cut -d ' ' -f2)
    motif_name=$(echo ${p} | cut -d ' ' -f3)
    fimo --o ../output/SPCG_promoters/fly_factor_survey/${motif_id}---${motif_name} --thresh 0.001 --motif ${motif_id} ../input/motif_databases.12.23/FLY/fly_factor_survey.meme ../output/SPCG_promoters/spcg_promoter_seq.fasta
done < ../output/SPCG_promoters/fly_factor_survey/motif_database.txt