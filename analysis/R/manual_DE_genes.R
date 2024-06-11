library(Seurat)
library(stringr)
library(singleCellNet)
library(Matrix)
library(pheatmap)
library(enrichR)
library(RColorBrewer)

enrichR::setEnrichrSite("FlyEnrichr")

TARGET_dir = file.path("results", ANALYSIS_VERSION, "manual_DE_genes")
dir.create(TARGET_dir)

wt_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13", "manual_celltype_object2.rds"))
rib_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_rib", "manual_celltyping_4_obj.rds"))

common_celltypes = intersect(wt_object@meta.data$manual_celltypes, rib_object@meta.data$manual_celltypes)
common_celltypes = common_celltypes[common_celltypes != 'Unknown']

# 
withr::with_dir(TARGET_dir, {
  for(celltype in common_celltypes) {
    print(celltype)
    dir.create(stringr::str_replace_all(celltype, "/", "-"))
    sub_wt_obj = subset(wt_object, subset = manual_celltypes == celltype) 
    sub_wt_obj@meta.data$experimental_condition = 'Wt'
    
    sub_rib_obj = subset(rib_object, subset = manual_celltypes == celltype) 
    sub_rib_obj@meta.data$experimental_condition = 'Rib'
    
    combined_obj = merge(sub_wt_obj, sub_rib_obj)
    object = Seurat::CreateSeuratObject(combined_obj@assays$RNA@counts, project = 'find_diff')
    object@meta.data$experimental_condition = combined_obj@meta.data$experimental_condition
    object@meta.data$batch = paste0(combined_obj@meta.data$experimental_condition, "_", combined_obj@meta.data$batch)
    
    object = Seurat::NormalizeData(object)
    markers <- FindMarkers(object = object, ident.1 = "Wt", logfc.threshold = 0, min.pct = 0, group.by = 'experimental_condition', test.use = 'wilcox') 
    write.csv(markers, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'DE_genes_wt_wilcox.csv'))
    
    markers = markers[markers$p_val_adj < 0.05, ]
    wt_markers = markers[markers$avg_log2FC > 0, ]
    rib_markers = markers[markers$avg_log2FC < 0, ]
    
    write.csv(wt_markers, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'markers_wt.csv'))
    write.csv(rib_markers, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'markers_rib.csv'))
    
    enrichment_results = enrichR::enrichr(
      genes = rownames(wt_markers), 
      databases = c(
        "GO_Biological_Process_2018", 
        "GO_Molecular_Function_2018", 
        'KEGG_2019'
      )
    )
    
    biological_analysis = enrichment_results$GO_Biological_Process_2018
    molecular_analysis = enrichment_results$GO_Molecular_Function_2018
    
    biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
    
    write.csv(biological_analysis, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'wt_biological_GO.csv'))
    write.csv(molecular_analysis, file = file.path(stringr::str_replace_all(celltype, "/", "-"), "wt_molecular_function_GO.csv"))
    
    # this is to get the rib 
    enrichment_results = enrichR::enrichr(
      genes = rownames(rib_markers), 
      databases = c(
        "GO_Biological_Process_2018", 
        "GO_Molecular_Function_2018", 
        'KEGG_2019'
      )
    )
    
    biological_analysis = enrichment_results$GO_Biological_Process_2018
    molecular_analysis = enrichment_results$GO_Molecular_Function_2018
    
    biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
    biological_analysis = biological_analysis[order(biological_analysis$Overlap, decreasing = TRUE), ]
    write.csv(biological_analysis, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'rib_biological_GO.csv'))
    write.csv(molecular_analysis, file = file.path(stringr::str_replace_all(celltype, "/", "-"), "rib_molecular_function_GO.csv"))
    
    ##############
    # might as well compute them in UMAP 
    object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
    object %<>% ScaleData(features = VariableFeatures(object = object))
    object %<>% RunPCA(features = VariableFeatures(object = object), npcs = 100)
    pc_cutoff = RMTstat::qmp( 1, ndf=length(Cells(object)), pdim=length(VariableFeatures(object)), var=1)
    singular_values = slot(Reductions(object, slot = "pca"), "stdev")
    is_significant = singular_values^2 > pc_cutoff
    num_pc = sum(is_significant)
    
    object %<>% RunPCA(features = VariableFeatures(object = object), npcs = num_pc)
    object %<>% FindNeighbors(dims = 1:num_pc)
    object %<>% RunUMAP(dim = 1:num_pc)
    object %<>% FindClusters()
    
    DimPlot(object)
    ggsave(filename = file.path(stringr::str_replace_all(celltype, "/", "-"), 'naive_cluster.png'), width = 8, height = 6)
    
    DimPlot(object, group.by = 'experimental_condition')
    ggsave(filename = file.path(stringr::str_replace_all(celltype, "/", "-"), 'naive_experimental_condition.png'), width = 8, height = 6)
    
    DimPlot(object, group.by = 'batch')
    ggsave(filename = file.path(stringr::str_replace_all(celltype, "/", "-"), 'naive_batch.png'), width = 8, height = 6)
    saveRDS(object, file = file.path(stringr::str_replace_all(celltype, "/", "-"), "ct_object.rds"))
  }
})

# reanalyze the enrichment with higher log2FC
withr::with_dir(TARGET_dir, {
  for(celltype in common_celltypes) {
    print(celltype)
    markers = read.csv(file.path(stringr::str_replace_all(celltype, "/", "-"), 'DE_genes_wt_wilcox.csv'), row.names = 1)
    
    markers = markers[markers$p_val_adj < 0.05, ]
    wt_markers = markers[markers$avg_log2FC > 0.25, ]
    rib_markers = markers[markers$avg_log2FC < -0.25, ]
    
    write.csv(wt_markers, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'markers_wt.csv'))
    write.csv(rib_markers, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'markers_rib.csv'))
    
    enrichment_results = enrichR::enrichr(
      genes = rownames(wt_markers), 
      databases = c(
        "GO_Biological_Process_2018", 
        "GO_Molecular_Function_2018", 
        'KEGG_2019'
      )
    )
    
    biological_analysis = enrichment_results$GO_Biological_Process_2018
    molecular_analysis = enrichment_results$GO_Molecular_Function_2018
    
    biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
    
    write.csv(biological_analysis, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'wt_biological_GO.csv'))
    write.csv(molecular_analysis, file = file.path(stringr::str_replace_all(celltype, "/", "-"), "wt_molecular_function_GO.csv"))
    
    # this is to get the rib 
    enrichment_results = enrichR::enrichr(
      genes = rownames(rib_markers), 
      databases = c(
        "GO_Biological_Process_2018", 
        "GO_Molecular_Function_2018", 
        'KEGG_2019'
      )
    )
    
    biological_analysis = enrichment_results$GO_Biological_Process_2018
    molecular_analysis = enrichment_results$GO_Molecular_Function_2018
    
    biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
    biological_analysis = biological_analysis[order(biological_analysis$Overlap, decreasing = TRUE), ]
    write.csv(biological_analysis, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'rib_biological_GO.csv'))
    write.csv(molecular_analysis, file = file.path(stringr::str_replace_all(celltype, "/", "-"), "rib_molecular_function_GO.csv"))
    
  }
})


