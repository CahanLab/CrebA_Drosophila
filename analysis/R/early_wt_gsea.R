TARGET_dir = file.path("results", ANALYSIS_VERSION, "early_wt_gsea")
dir.create(TARGET_dir)

object = readRDS(file.path("results", ANALYSIS_VERSION, 'harmonized_wildtype_data/stage10-12_reharmonized_seurat.rds'))
pathway_list = readRDS('accessory_data/GO_Biological_Processes_2018/GO_Biological_Process.rds')

##### compute GSEA #####
withr::with_dir(TARGET_dir, {
  for(celltype in unique(object$manual_celltypes)) {
    if(stringr::str_replace_all(celltype, "/", "-") %in% list.dirs(".", recursive = FALSE, full.names = FALSE)) {
      next() # if we have already calculated the marker genes, then we skip 
    } else {
      print(stringr::str_replace_all(celltype, "/", "-"))
      dir.create(stringr::str_replace_all(celltype, "/", "-"))
      markers = SeuratWrappers::RunPresto(object, ident.1 = celltype, logfc.threshold = 0, min.pct = 0.1, group.by = 'manual_celltypes') 
      write.csv(markers, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'markers_genes.csv'))
      
      ranks <- markers$avg_log2FC
      names(ranks) <- rownames(markers)
      
      # allPathways is predefined list of signatures
      fgseaRes <- fgsea(pathways = pathway_list, 
                        stats = ranks,
                        minSize=10,
                        maxSize=500)
      fgseaRes = data.frame(fgseaRes)
      fgseaRes <- apply(fgseaRes,2,as.character)
      
      write.csv(fgseaRes, file = file.path(file.path(stringr::str_replace_all(celltype, "/", "-")), 'gsea_results_early_wt.csv'))
      
      markers = markers[markers$p_val_adj < 0.05, ]
      markers = markers[markers$avg_log2FC > 0, ]
      
      enrichment_results = enrichR::enrichr(
        genes = rownames(markers), 
        databases = c(
          "GO_Biological_Process_2018", 
          "GO_Molecular_Function_2018"
        )
      )
      
      biological_analysis = enrichment_results$GO_Biological_Process_2018
      molecular_analysis = enrichment_results$GO_Molecular_Function_2018
      
      biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
      molecular_analysis = molecular_analysis[molecular_analysis$Adjusted.P.value < 0.05, ]
      write.csv(biological_analysis, file = file.path(stringr::str_replace_all(celltype, "/", "-"), 'biological_GO.csv'))
      write.csv(molecular_analysis, file = file.path(stringr::str_replace_all(celltype, "/", "-"), "molecular_function_GO.csv"))
    }
  }
})

##### make preliminary plots #####
all_folders = list.dirs(TARGET_dir, recursive = FALSE)
for(my_folder in all_folders) {
  print(my_folder)
  gsea_results = read.csv(file.path(my_folder, 'gsea_results_early_wt.csv'), row.names = 1)
  gsea_results = gsea_results[gsea_results$padj < 0.05, ]
  gsea_results = gsea_results[gsea_results$NES > 0, ]
  gsea_results$num_leadingEdges = NULL
  for(temp_index in rownames(gsea_results)) {
    leading_edge = eval(parse(text = gsea_results[temp_index, 'leadingEdge']))
    gsea_results[temp_index, 'num_leadingEdges'] = length(leading_edge)
  }
  gsea_results$GeneRatio = gsea_results$num_leadingEdges / gsea_results$size
  gsea_results$NES = as.numeric(gsea_results$NES)
  gsea_results = gsea_results[order(gsea_results$NES, decreasing = TRUE), ]
  min_index = min(nrow(gsea_results), 20)
  gsea_results = gsea_results[1:min_index, ]
  
  p = ggplot(gsea_results, aes(x = reorder(pathway, -padj), y = NES)) + 
    geom_point(aes(size = GeneRatio, color = padj)) +
    scale_size_continuous(range = c(4,8)) +
    theme_bw(base_size = 14) +
    scale_colour_gradient(limits=c(0, 0.05), low="red") +
    ylab('Normalized Enrichment Scores') +
    xlab("GO terms") +
    ggtitle(stringr::str_split(my_folder, "/")[[1]][length(stringr::str_split(my_folder, "/")[[1]])]) + 
    coord_flip()
  ggsave(filename = file.path(my_folder, 'GO_enrichment_dot.png'), plot = p, width = 12, height = 10)
}
