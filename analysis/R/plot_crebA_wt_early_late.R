TARGET_dir = file.path("results", ANALYSIS_VERSION, "plot_crebA_wt_early_late")
dir.create(TARGET_dir, recursive = TRUE)

wt_late_object = readRDS("accessory_data/wild_type_seurats/stage_13-16/manual_celltype_object4.rds")
wt_early_object = readRDS("accessory_data/wild_type_seurats/stage_10-12/manual_celltype_object1.rds")

crebA_late_object = readRDS("results/v19/manual_annotation_crebA/manual_celltype_object.rds")
crebA_early_object = readRDS("results/v19/manual_annotation_crebA_early/manual_celltype_object.rds")

common_cts_list = Reduce(intersect, list(unique(wt_early_object@meta.data$manual_celltypes),
                                    unique(wt_late_object@meta.data$manual_celltypes), 
                                    unique(crebA_early_object@meta.data$manual_celltypes), 
                                    unique(crebA_late_object@meta.data$manual_celltypes)))
for(common_ct in common_cts_list) {
  print(common_ct) 
  sub_wt_late_object = subset(wt_late_object, subset = manual_celltypes == common_ct) 
  sub_wt_early_object = subset(wt_early_object, subset = manual_celltypes == common_ct) 
  sub_crebA_late_object = subset(crebA_late_object, subset = manual_celltypes == common_ct) 
  sub_crebA_early_object = subset(crebA_early_object, subset = manual_celltypes == common_ct) 
  
  sub_wt_late_object@meta.data$experimental_condition = 'wt_late'
  sub_wt_early_object@meta.data$experimental_condition = 'wt_early'
  sub_crebA_late_object@meta.data$experimental_condition = 'crebA_late'
  sub_crebA_early_object@meta.data$experimental_condition = 'crebA_early'
  
  combined_obj = merge(x = sub_wt_late_object, y = c(sub_wt_early_object, sub_crebA_late_object, sub_crebA_early_object))
  object = Seurat::CreateSeuratObject(combined_obj@assays$RNA@counts, project = 'combine')
  object@meta.data$experimental_condition = combined_obj@meta.data$experimental_condition
  object@meta.data$batch = paste0(combined_obj@meta.data$experimental_condition, "_", combined_obj@meta.data$batch)
  
  object@meta.data$batch_condition = paste0(object@meta.data$batch, "-", object@meta.data$experimental_condition)  
  
  object = Seurat::NormalizeData(object)
  object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  object %<>% ScaleData(features = VariableFeatures(object = object))
  object %<>% RunPCA(features = VariableFeatures(object = object), npcs = 100)
  
  num_pc = 40 # this is hard coded 
  # after getting the significant PCs, rerun PCA 
  object %<>% RunPCA(features = VariableFeatures(object = object), npcs = num_pc)
  object %<>% FindNeighbors(dims = 1:num_pc, k.param = 10) # to really break down the clusters. Probably will need to do the same for the wild types
  object %<>% RunUMAP(dim = 1:num_pc)
  
  dir.create(file.path(TARGET_dir, common_ct), recursive = TRUE)
  p = DimPlot(object, group.by = "experimental_condition", label = TRUE, label.size = 4)
  ggsave(filename = file.path(TARGET_dir, common_ct, 'ct_UMAP.png'), plot = p, width = 7, height = 6)
  saveRDS(object, file = file.path(TARGET_dir, common_ct, 'subset_object.rds'))
}
