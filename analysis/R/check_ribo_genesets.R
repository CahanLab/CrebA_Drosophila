library(Seurat)
library(stringr)
library(singleCellNet)
library(Matrix)
library(pheatmap)
library(enrichR)
library(RColorBrewer)

enrichR::setEnrichrSite("FlyEnrichr")

TARGET_dir = file.path("results", ANALYSIS_VERSION, "check_ribo_genesets")
dir.create(TARGET_dir)

wt_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13", "manual_celltype_object2.rds"))
rib_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_rib", "manual_celltyping_4_obj.rds"))

common_celltypes = intersect(wt_object@meta.data$manual_celltypes, rib_object@meta.data$manual_celltypes)
common_celltypes = common_celltypes[common_celltypes != 'Unknown']

ribo_genes_df = readxl::read_xlsx("accessory_data/genesets/ribo_genes/ribo_genes.xlsx")

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
    
    all_ribos_fly = intersect(ribo_genes_df$`FlyBase ID`, names(all_genes))
    all_ribos_genes = as.vector(all_genes[all_ribos_fly])
    object <- AddModuleScore(object,
                           features = list(all_ribos_genes),
                           name="all_ribo_genes")
    
    VlnPlot(object, features = paste0('all_ribo_genes', 1), group.by = 'experimental_condition', pt.size = 0)
    ggsave(filename = file.path(stringr::str_replace_all(celltype, "/", "-"), 'ribo_gene_exp_condition.png'), width = 8, height = 6)
    
    VlnPlot(object, features = paste0('all_ribo_genes', 1), group.by = 'batch', pt.size = 0)
    ggsave(filename = file.path(stringr::str_replace_all(celltype, "/", "-"), 'ribo_gene_exp_batch.png'), width = 8, height = 6)
    
    SG_ribos_fly = intersect(as.vector(ribo_genes_df$`FlyBase ID`[ribo_genes_df$`Bound by Rib (SG)` == 'Yes']), names(all_genes))
    SG_ribos_genes = as.vector(all_genes[SG_ribos_fly])
    object <- AddModuleScore(object,
                             features = list(SG_ribos_genes),
                             name="SG_ribo_genes")
    
    VlnPlot(object, features = paste0('SG_ribo_genes', 1), group.by = 'experimental_condition', pt.size = 0)
    ggsave(filename = file.path(stringr::str_replace_all(celltype, "/", "-"), 'SG_ribo_gene_exp_condition.png'), width = 8, height = 6)
    
    VlnPlot(object, features = paste0('SG_ribo_genes', 1), group.by = 'batch', pt.size = 0)
    ggsave(filename = file.path(stringr::str_replace_all(celltype, "/", "-"), 'SG_ribo_gene_exp_batch.png'), width = 8, height = 6)
    
    TS_ribos_fly = intersect(as.vector(ribo_genes_df$`FlyBase ID`[ribo_genes_df$`Bound by Rib (Tr.)` == 'Yes']), names(all_genes))
    TS_ribos_genes = as.vector(all_genes[TS_ribos_fly])
    object <- AddModuleScore(object,
                             features = list(TS_ribos_genes),
                             name="TS_ribo_genes")
    
    VlnPlot(object, features = paste0('TS_ribo_genes', 1), group.by = 'experimental_condition', pt.size = 0)
    ggsave(filename = file.path(stringr::str_replace_all(celltype, "/", "-"), 'TS_ribo_gene_exp_condition.png'), width = 8, height = 6)
    
    VlnPlot(object, features = paste0('TS_ribo_genes', 1), group.by = 'batch', pt.size = 0)
    ggsave(filename = file.path(stringr::str_replace_all(celltype, "/", "-"), 'TS_ribo_gene_exp_batch.png'), width = 8, height = 6)
    
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

