TARGET_dir = file.path("results", ANALYSIS_VERSION, "DE_genes_early_rib_wt")
dir.create(TARGET_dir)

wt_object = readRDS("accessory_data/wild_type_seurats/stage_10-12/manual_celltype_object1.rds")
mut_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_rib_early/manual_celltype_object.rds"))

common_celltypes = intersect(wt_object@meta.data$manual_celltypes, mut_object@meta.data$manual_celltypes)
pathway_list = readRDS('accessory_data/GO_Biological_Processes_2018/GO_Biological_Process.rds')

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
    rank_sum_test = presto::wilcoxauc(object, group_by = "experimental_condition")
    rank_sum_test = rank_sum_test[rank_sum_test$pct_in > 10 | rank_sum_test$pct_out > 10, ]
    
    # compute the DE genes for wild type 
    sub_rank_sum_test = rank_sum_test[rank_sum_test$group == 'Wt', ]
    write.csv(sub_rank_sum_test, file = file.path(stringr::str_replace_all(celltype, "/", "-"), "wt_DE_genes.csv"))
    
    ranks <- sub_rank_sum_test$logFC
    names(ranks) <- sub_rank_sum_test$feature
    fgseaRes <- fgsea::fgsea(pathways = pathway_list, 
                      stats = ranks,
                      minSize=10,
                      maxSize=500)
    
    fgseaRes = data.frame(fgseaRes) 
    fgseaRes = apply(fgseaRes,2,as.character)
    fgseaRes = as.data.frame(fgseaRes)
    fgseaRes = fgseaRes[!is.na(fgseaRes$padj), ]
    fgseaRes = fgseaRes[!is.na(fgseaRes$NES), ]
    
    write.csv(fgseaRes, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'wt_gsea_results.csv'))
    
    # this is to calculate the mutants 
    sub_rank_sum_test = rank_sum_test[rank_sum_test$group == 'mut', ]
    write.csv(sub_rank_sum_test, file = file.path(stringr::str_replace_all(celltype, "/", "-"), "mut_DE_genes.csv"))
    
    ranks <- sub_rank_sum_test$logFC
    names(ranks) <- sub_rank_sum_test$feature
    fgseaRes <- fgsea::fgsea(pathways = pathway_list, 
                             stats = ranks,
                             minSize=10,
                             maxSize=500)
    
    fgseaRes = data.frame(fgseaRes) 
    fgseaRes = apply(fgseaRes,2,as.character)
    fgseaRes = as.data.frame(fgseaRes)
    fgseaRes = fgseaRes[!is.na(fgseaRes$padj), ]
    fgseaRes = fgseaRes[!is.na(fgseaRes$NES), ]
    
    write.csv(fgseaRes, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'mut_gsea_results.csv'))
    
  }
})
