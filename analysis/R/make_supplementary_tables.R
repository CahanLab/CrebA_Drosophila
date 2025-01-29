TARGET_dir = file.path("results", ANALYSIS_VERSION, "SuppTabs")

##### make the cluster DE genes #####
DE_tab = read.csv(file.path("results", ANALYSIS_VERSION, "crebA_early_integrated/marker_genes.csv"), row.names = 1)
write.csv(DE_tab, file.path(TARGET_dir, 'early_CrebA_marker_genes.csv'))

# most likely not going to be used in the main analysis 
# but just generate this in case we need to
DE_tab = read.csv(file.path('results', ANALYSIS_VERSION, 'crebA_rep3/marker_genes.csv'), row.names = 1)
write.csv(DE_tab, file.path(TARGET_dir, 'CrebA_marker_genes.csv'))

##### make the DE genes between wt and CrebA ###### 
input_path = file.path("results", ANALYSIS_VERSION, "DE_genes_early_crebA_wt/")

excel_list = list()
for(tmp_ct in list.files(input_path)) {
  tmp_DE = read.csv(file.path(input_path, tmp_ct, "mut_DE_genes.csv"), row.names = 1)
  tmp_DE = tmp_DE[order(tmp_DE$logFC), ]
  excel_list[[tmp_ct]] = tmp_DE
}
openxlsx::write.xlsx(excel_list, file = file.path(TARGET_dir, "DE_genes_early_CrebA_wt.xlsx"))

# this is for the late CrebA 
input_path = file.path("results", ANALYSIS_VERSION, "DE_genes_crebA_wt/")

excel_list = list()
for(tmp_ct in list.files(input_path)) {
  tmp_DE = read.csv(file.path(input_path, tmp_ct, "mut_DE_genes.csv"), row.names = 1)
  tmp_DE = tmp_DE[order(tmp_DE$logFC), ]
  excel_list[[tmp_ct]] = tmp_DE
}
openxlsx::write.xlsx(excel_list, file = file.path(TARGET_dir, "DE_genes_CrebA_wt.xlsx"))

##### make the GSEA tables between early WT and CrebA #####
input_path = file.path("results", ANALYSIS_VERSION, "DE_genes_early_crebA_wt/")

excel_list = list()
for(tmp_ct in list.files(input_path)) {
  tmp_df = read.csv(file.path(input_path, tmp_ct, "mut_gsea_results.csv"), row.names = 1)
  tmp_df = tmp_df[tmp_df$NES > 0, ]
  tmp_df = tmp_df[order(tmp_df$padj), ]
  excel_list[[tmp_ct]] = tmp_DE
}
openxlsx::write.xlsx(excel_list, file = file.path(TARGET_dir, "mut_GSEA_early_CrebA_wt.xlsx"))

# get the wildtype enrichment 
input_path = file.path("results", ANALYSIS_VERSION, "DE_genes_early_crebA_wt/")

excel_list = list()
for(tmp_ct in list.files(input_path)) {
  tmp_df = read.csv(file.path(input_path, tmp_ct, "wt_gsea_results.csv"), row.names = 1)
  tmp_df = tmp_df[tmp_df$NES > 0, ]
  tmp_df = tmp_df[order(tmp_df$padj), ]
  excel_list[[tmp_ct]] = tmp_DE
}
openxlsx::write.xlsx(excel_list, file = file.path(TARGET_dir, "wt_GSEA_early_CrebA_wt.xlsx"))
