#/bin/bash 
mkdir ../output/SPCG_promoters/CrebA_motifs/
cat ../input/previous_CrebA_motifs/CrebAMotifs_meme.meme | grep MOTIF > ../output/SPCG_promoters/CrebA_motifs/motif_database.txt 
while read p; do
    motif_id=$(echo ${p} | cut -d ' ' -f2)
    motif_name=$(echo ${p} | cut -d ' ' -f3)
    echo ${motif_name}
    fimo --o ../output/SPCG_promoters/CrebA_motifs/${motif_name} --thresh 0.001 --motif ${motif_id} ../input/previous_CrebA_motifs/CrebAMotifs_meme.meme ../output/SPCG_promoters/spcg_promoter_seq.fasta
done < ../output/SPCG_promoters/CrebA_motifs/motif_database.txt 

