TARGET_dir = file.path("results", ANALYSIS_VERSION, "early_crebA_wt_ct_objects")
dir.create(TARGET_dir)

wt_object = readRDS(file.path('results', ANALYSIS_VERSION, "harmonized_wildtype_data/stage10-12_reharmonized_seurat.rds"))
wt_object@meta.data$manual_celltypes = wt_object@meta.data$new_celltypes

mut_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_crebA_early/manual_celltype_object.rds"))

common_celltypes = intersect(wt_object@meta.data$manual_celltypes, mut_object@meta.data$manual_celltypes)

withr::with_dir(TARGET_dir, {
  for(celltype in common_celltypes) {
    dir.create(stringr::str_replace_all(celltype, "/", "-"))
    sub_wt_obj = subset(wt_object, subset = manual_celltypes == celltype) 
    sub_wt_obj@meta.data$experimental_condition = 'Wt'
    
    sub_mut_obj = subset(mut_object, subset = manual_celltypes == celltype) 
    sub_mut_obj@meta.data$experimental_condition = 'mut'
    
    combined_obj = merge(sub_wt_obj, sub_mut_obj)
    object = Seurat::CreateSeuratObject(combined_obj@assays$RNA@counts, project = 'find_diff')
    object@meta.data$experimental_condition = combined_obj@meta.data$experimental_condition
    object@meta.data$batch = paste0(combined_obj@meta.data$experimental_condition, "_", combined_obj@meta.data$batch)
    
    object = Seurat::NormalizeData(object)
    
    object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
    object %<>% ScaleData(features = VariableFeatures(object = object))
    object %<>% RunPCA(features = VariableFeatures(object = object), npcs = 100)
    
    num_pc = 60 # this is hard coded 
    object %<>% RunPCA(features = VariableFeatures(object = object), npcs = num_pc)
    object %<>% FindNeighbors(dims = 1:num_pc, k.param = 10) # to really break down the clusters. Probably will need to do the same for the wild types
    object %<>% RunUMAP(dim = 1:num_pc)
    object %<>% FindClusters()
    
    p = DimPlot(object, group.by = "batch")
    ggsave(filename = file.path(stringr::str_replace_all(celltype, "/", "-"), 'batch.png'), plot = p, width = 8, height = 8)
    p = DimPlot(object, group.by = 'experimental_condition')
    ggsave(filename = file.path(stringr::str_replace_all(celltype, "/", "-"), 'experimental_conditions.png'), plot = p, width = 8, height = 8)
    
    object %<>% harmony::RunHarmony("batch")
    object %<>% FindNeighbors(reduction = "harmony", dims = 1:num_pc, k.param = 10)
    object %<>% RunUMAP(reduction = "harmony", dim = 1:num_pc)
    
    p = DimPlot(object, group.by = "batch")
    ggsave(filename = file.path(stringr::str_replace_all(celltype, "/", "-"), 'harmony_batch.png'), plot = p, width = 8, height = 8)
    p = DimPlot(object, group.by = 'experimental_condition')
    ggsave(filename = file.path(stringr::str_replace_all(celltype, "/", "-"), 'harmony_experimental_conditions.png'), plot = p, width = 8, height = 8)
    
    saveRDS(object, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'object.rds'))
    
  }
})
