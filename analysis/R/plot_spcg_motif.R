TARGET_dir = file.path("results", ANALYSIS_VERSION, "Figures/plot_spcg_motif")
dir.create(TARGET_dir)

##### look at SPCGs and CrebA motif #####
creba_motif = read.csv("accessory_data/CrebA_ChIP_Analysis/output/SPCG_promoters/crebA_hits_0.0001.csv", row.names = 1)
creba_sum = apply(creba_motif, MARGIN = 1, FUN = sum) / 90
sum(apply(creba_motif, MARGIN = 2, FUN = sum) == 0) / 90

# seems like the number is a bit low. Let me rethink about the thresholding again 
# might have to adjust the p value to 0.001 --> what was done before 

