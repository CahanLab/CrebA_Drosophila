#!/bin/bash

# get the early data 
Rscript --vanilla ../accessory_data/FlyPhoneDB/FlyPhone_parallel_batch.R \
		--matrix ../results/v19/FlyPhoneDB_run/early_wt_count_matrix.csv \
		--metadata ../results/v19/FlyPhoneDB_run/early_wt_sampTab.csv \
		--lrpair ../accessory_data/FlyPhoneDB/annotation/Ligand_receptor_pair_high_confident_2021vs1_clean.txt \
		--corecomponents ../accessory_data/FlyPhoneDB/annotation/Pathway_core_components_2021vs1_clean.txt \
		--cores 8 \
		--output ../results/v19/FlyPhoneDB_run/early_wt/

# get the late data 
Rscript --vanilla ../accessory_data/FlyPhoneDB/FlyPhone_parallel_batch.R \
		--matrix ../results/v19/FlyPhoneDB_run/late_wt_count_matrix.csv \
		--metadata ../results/v19/FlyPhoneDB_run/late_wt_sampTab.csv \
		--lrpair ../accessory_data/FlyPhoneDB/annotation/Ligand_receptor_pair_high_confident_2021vs1_clean.txt \
		--corecomponents ../accessory_data/FlyPhoneDB/annotation/Pathway_core_components_2021vs1_clean.txt \
		--cores 8 \
		--output ../results/v19/FlyPhoneDB_run/late_wt/

# get the early crebA data 
Rscript --vanilla ../accessory_data/FlyPhoneDB/FlyPhone_parallel_batch.R \
		--matrix ../results/v19/FlyPhoneDB_run/early_crebA_count_matrix.csv \
		--metadata ../results/v19/FlyPhoneDB_run/early_crebA_sampTab.csv \
		--lrpair ../accessory_data/FlyPhoneDB/annotation/Ligand_receptor_pair_high_confident_2021vs1_clean.txt \
		--corecomponents ../accessory_data/FlyPhoneDB/annotation/Pathway_core_components_2021vs1_clean.txt \
		--cores 8 \
		--output ../results/v19/FlyPhoneDB_run/early_crebA/