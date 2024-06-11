library(fgsea)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "manual_DE_genes")
all_paths = list.dirs(TARGET_dir, recursive = FALSE)

pathway_list = readRDS('accessory_data/GO_Biological_Processes_2018/GO_Biological_Process.rds')
for(current_path in all_paths) {
  print(current_path)
  DE_genes = read.csv(file.path(current_path, 'DE_genes_wt_wilcox.csv'), row.names = 1)
  sub_DE_genes = DE_genes[DE_genes$pct.1 > 0.1, ]
  ranks <- sub_DE_genes$avg_log2FC
  names(ranks) <- rownames(sub_DE_genes)
  
  # allPathways is predefined list of signatures
  fgseaRes <- fgsea(pathways = pathway_list, 
                    stats = ranks,
                    minSize=10,
                    maxSize=500,
                    nperm=1000000)
  fgseaRes = data.frame(fgseaRes)
  fgseaRes <- apply(fgseaRes,2,as.character)
  
  write.csv(fgseaRes, file = file.path(current_path, 'gsea_results_wt.csv'))
  
  # allPathways is predefined list of signatures
  sub_DE_genes = DE_genes[DE_genes$pct.2 > 0.1, ]
  ranks <- sub_DE_genes$avg_log2FC
  names(ranks) <- rownames(sub_DE_genes)
  fgseaRes <- fgsea(pathways = pathway_list, 
                    stats = -ranks,
                    minSize=10,
                    maxSize=500,
                    nperm=1000000)
  fgseaRes = data.frame(fgseaRes)
  fgseaRes <- apply(fgseaRes,2,as.character)
  
  write.csv(fgseaRes, file = file.path(current_path, 'gsea_results_rib.csv'))
}

all_folders = list.dirs(TARGET_dir, recursive = FALSE)
for(my_folder in all_folders) {
  print(my_folder)
  gsea_results = read.csv(file.path(my_folder, 'gsea_results_rib.csv'), row.names = 1)
  gsea_results = gsea_results[gsea_results$padj < 0.05, ]
  gsea_results$num_leadingEdges = NULL
  for(temp_index in rownames(gsea_results)) {
    leading_edge = eval(parse(text = gsea_results[temp_index, 'leadingEdge']))
    gsea_results[temp_index, 'num_leadingEdges'] = length(leading_edge)
  }
  gsea_results$GeneRatio = gsea_results$num_leadingEdges / gsea_results$size
  gsea_results$NES = as.numeric(gsea_results$NES)
  gsea_results = gsea_results[order(abs(gsea_results$NES), decreasing = TRUE), ]
  min_index = min(nrow(gsea_results), 20)
  gsea_results = gsea_results[1:min_index, ]
  
  p = ggplot(gsea_results, aes(x = reorder(pathway, -padj), y = NES)) + 
    geom_point(aes(size = GeneRatio, color = padj)) +
    scale_size_continuous(range = c(4,8)) +
    theme_bw(base_size = 14) +
    scale_colour_gradient(limits=c(0, 0.05), low="red") +
    ylab('Normalized Enrichment Scores') +
    xlab("GO terms") +
    ggtitle(paste0(stringr::str_split(my_folder, "/")[[1]][length(stringr::str_split(my_folder, "/")[[1]])], " effect of Rib knockout")) + 
    coord_flip()
  ggsave(filename = file.path(my_folder, 'GO_enrichment_dot.png'), plot = p, width = 12, height = 10)
}


