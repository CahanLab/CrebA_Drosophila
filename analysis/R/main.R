# To reproduce our work, run this script. 
# We used R version 4.1.2.

# Make it work on each of our computers
try({setwd("~/Dropbox/salivary_gland_ribbon/analysis")}, silent = T)
try({setwd("~/Dropbox (CahanLab)/salivary_gland_ribbon/analysis")}, silent = T)
try({setwd("c:/Users/DAndrewLab/Desktop/wt_scRNA/analysis/")}, silent = T)
renv::activate() 

ANALYSIS_VERSION = "v19" # Where to put output
metadata = read.table(header = T, text=
  "sample cellranger
  rib_rep1 2022-04-11_scRNA_10x_3prime
  rib_rep2 2022-06-01_scRNA_10x_3prime
   crebA_rep1 2022-12-07_scRNA_10x_3prime_CrebA2
   crebA_early_rep1 2022-12-07_scRNA_10x_3prime_CrebA
   rib_early_rep1 2022-12-13_scRNA_10x_3prime_Rib_early1
   crebA_early_rep2 2022-12-14_scRNA_10x_3prime_CrebA_early2
   crebA_rep2 2022-12-15_scRNA_10x_3prime_CrebA_late2
   rib_early_rep2 2022-12-15_scRNA_10x_3prime_Rib_early2
   crebA_rep3 2022-12-16_scRNA_10x_3prime_CrebA_late3
   crebA_rep4 2023-04-05_scRNA_10x_3prime_CrebA_late4
   crebA_rep5 2023-06-20_scRNA_10x_3prime_CrebA_late5")

SAMPLE = "all"
source("R/set_up_environment.R")

# change this to 1:14 when you are done with everything
# when i == 12 CrebA_rep2 is very weird. Look at it individually 
for(i in 1:11){ 
  SAMPLE = metadata$sample[i]  
  CELLRANGER = metadata$cellranger[i] 
  source("R/set_up_environment.R") # Uses SAMPLE to set location of output
  source("R/filtering_and_qc.R") 
  source("R/refined_clustering.R")
  source("R/annotate_clusters_automated_BDGP.R") # perform a quick annotation using genes from BDGP 
}

# make the folder path for the figures 
if(dir.exists(file.path("results", ANALYSIS_VERSION, 'Figures')) == FALSE) {
  dir.create(file.path("results", ANALYSIS_VERSION, 'Figures'))
}

# harmonize the wildtype cell types 
source("R/harmonize_wildtype.R")

# look at embryo level DE genes in terms of the embryo level 
# the purpose is to compare and contrast with microarray 
source("R/DE_genes_embryo_level.R")

# this is to look at the cell cell communication 
source("R/FlyPhoneDB_wt.R")

# integrate the early CrebA data 
source("R/integrate_early_crebA.R")

# manual label the early crebA data 
source("R/manual_label_early_crebA.R")

# find DE genes between different cell types in CrebA 
source("R/early_crebA_gsea.R")

# find the DE genes betweent different cell types in early wildtype 
source("R/early_wt_gsea.R")

# find DE genes between the same cell types 
source("R/DE_genes_early_crebA_wt.R")

# project same tissues onto UMAPs 
source("R/crebA_wt_early_ct_subset.R")

# we will just use CrebA rep 3 as the only CrebA sample. No replicate. 
source("R/manual_label_crebA.R")

# find DE genes between the same cell types 
source("R/DE_genes_crebA_wt.R")

##### make the plotting scripts ######

# make the color palette 
source("R/create_color_palette.R")

# make the plot for the cell types annotations in both CrebA and wildtype 
source("R/plot_celltypes.R")

# make the plots for the statics of CrebA dataset 
source("R/plot_crebA_stats.R")

# look at the spcg in the wild type data 
source("R/plot_spcg_wt.R")

# plot out the logFC of SPCGs among cell types 
source("R/plot_spcg_logFC.R")

# plot out the correlation between logFC and normalized expression of CrebA 
source("R/plot_logFC_norm_Exp_crebA.R")

# check interested genes 
source("R/check_genes_UMAP.R")

# look at the secretome in wildtype data 
source("R/plot_secretome_wt.R")

# plot some of the GSEA results 
source("R/plot_gsea_crebA_wt_early.R")

# plot some of the GSEA results 
source("R/plot_gsea_crebA_wt.R")

# plot out apodemes and optic lobe genes 
source("R/plot_ap_ol_genes.R")
###### unused scripts ######

# make the plots of the common cell types 
source("R/plot_crebA_wt_early_late.R")

# plot out the DE genes between different cell types in CrebA 
source("R/plot_early_crebA_gsea.R")
