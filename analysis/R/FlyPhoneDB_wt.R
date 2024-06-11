TARGET_dir = file.path("results", ANALYSIS_VERSION, "FlyPhoneDB_run")
dir.create(TARGET_dir)

##### curate the early seurat object #####
object = readRDS("accessory_data/wild_type_seurats/stage_10-12/manual_celltype_object1.rds")
counts_matrix = as.matrix(object@assays$RNA@counts)
write.csv(counts_matrix, file = file.path(TARGET_dir, "early_wt_count_matrix.csv"))
ct_tab = object@meta.data
new_df = subset(ct_tab, select = manual_celltypes)
colnames(new_df) = c("celltype")
write.csv(new_df, file = file.path(TARGET_dir, "early_wt_sampTab.csv"))

##### curate the late seurat object #####
object = readRDS("accessory_data/wild_type_seurats/stage_13-16/manual_celltype_object4.rds")
counts_matrix = as.matrix(object@assays$RNA@counts)
write.csv(counts_matrix, file = file.path(TARGET_dir, "late_wt_count_matrix.csv"))
ct_tab = object@meta.data
new_df = subset(ct_tab, select = manual_celltypes)
colnames(new_df) = c("celltype")
write.csv(new_df, file = file.path(TARGET_dir, "late_wt_sampTab.csv"))

##### curate the early creba seurat object #####
object = readRDS("results/v19/manual_annotation_crebA/manual_celltype_object.rds")
counts_matrix = as.matrix(object@assays$RNA@counts)
write.csv(counts_matrix, file = file.path(TARGET_dir, "early_crebA_count_matrix.csv"))
ct_tab = object@meta.data
new_df = subset(ct_tab, select = manual_celltypes)
colnames(new_df) = c("celltype")
write.csv(new_df, file = file.path(TARGET_dir, "early_crebA_sampTab.csv"))

###### contrast the difference between early wt and crebA #####
compiled_flyphoneDB <- function(flyphoneDB) { 
  non_int_cols = colnames(flyphoneDB)[(grepl("_pvalues", colnames(flyphoneDB)) | grepl("_score", colnames(flyphoneDB))) == FALSE]
  int_cols_pvals = colnames(flyphoneDB)[grepl('_pvalues', colnames(flyphoneDB))]
  int_cols_score = colnames(flyphoneDB)[grepl("_score", colnames(flyphoneDB))]
  
  big_df = data.frame()
  for(temp_pvals in int_cols_pvals) {
    temp_df = flyphoneDB[, non_int_cols]
    temp_interactions = stringr::str_remove(string = temp_pvals, pattern = '_pvalues')
    temp_score = int_cols_score[grepl(temp_interactions, int_cols_score)]
    temp_df$flyphone_scores = flyphoneDB[, temp_score]
    temp_df$flyphone_pvals = flyphoneDB[, temp_pvals]
    temp_df$ct_interactions = temp_interactions
    big_df = rbind(big_df, temp_df)
  }
  
  return(big_df)
}

wt_flyphone = read.csv("results/v19/FlyPhoneDB_run/early_wt/interaction_list.csv", row.names = 1)
wt_flyphone_long = compiled_flyphoneDB(wt_flyphone)
wt_flyphone_long = wt_flyphone_long[wt_flyphone_long$flyphone_pvals <= 0.05, ]
wt_flyphone_long$flyphone_scores = NULL
write.csv(wt_flyphone_long, file = file.path(TARGET_dir, "sig_early_wt_interactions.csv"))

crebA_flyphone = read.csv("results/v19/FlyPhoneDB_run/early_crebA/interaction_list.csv", row.names = 1)
crebA_flyphone_long = compiled_flyphoneDB(crebA_flyphone)
crebA_flyphone_long = crebA_flyphone_long[crebA_flyphone_long$flyphone_pvals <= 0.05, ]
crebA_flyphone_long$flyphone_scores = NULL
write.csv(crebA_flyphone_long, file = file.path(TARGET_dir, "sig_early_crebA_interactions.csv"))





