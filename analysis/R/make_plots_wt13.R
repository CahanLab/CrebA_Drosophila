library(ggplot2)
library(RColorBrewer)
library(ggdendroplot)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "wt13_plots")
dir.create(TARGET_dir)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
qual_col_pals = qual_col_pals[rownames(qual_col_pals) != 'Dark2', ]
qual_col_pals = qual_col_pals[rownames(qual_col_pals) != 'Set2', ]

col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = col_vector[col_vector != col_vector[4]]
col_vector = col_vector[col_vector != col_vector[4]]

wt_object = readRDS(file.path('results', ANALYSIS_VERSION, 'manual_annotation_wt13/manual_celltype_object3.rds'))

# plot out the cell type UMAP 
p = DimPlot(wt_object, group.by = 'manual_celltypes', label = TRUE, label.size = 5) +
  ggtitle("Celltype Labels") + 
  xlim(c(-8, 12))
ggsave(file.path(TARGET_dir, "cell_type_UMAP/UMAP.png"), plot = p, width = 20, height = 10)

# plot out the batch to demonstrate the batch effect 
unique(wt_object@meta.data$batch)
wt_object@meta.data$batch_labels = NA
wt_object@meta.data[wt_object@meta.data$batch == 'rep_1', 'batch_labels'] = 'batch_1'
wt_object@meta.data[wt_object@meta.data$batch == 'rep_3', 'batch_labels'] = 'batch_2'

p = DimPlot(wt_object, group.by = 'batch_labels', label.size = 5) +
  ggtitle("Batch Label") + 
  xlim(c(-8, 12))
ggsave(file.path(TARGET_dir, "batch_label_UMAP/UMAP.png"), plot = p, width = 12, height = 8)


# plot out the cell clusters 
p = DimPlot(wt_object, group.by = 'seurat_clusters', label = TRUE, label.size = 5) +
  ggtitle("Cluster Labels")
ggsave(file.path(TARGET_dir, "cluster_UMAP/UMAP.png"), plot = p, width = 12, height = 8)

# plot out the 3D UMAP 
wt_object %<>% RunUMAP(reduction = "harmony", dim = 1:66, n.components = 3L)
wt_object[["umap3D"]] = Reductions(wt_object, "umap")

dir.create(file.path(TARGET_dir, "3D_UMAP"))
withr::with_dir(file.path(TARGET_dir, "3D_UMAP"), {
  save3DUMAPCategorical(wt_object, group.by = 'manual_celltypes', umapkey = "umap3D")
})

# plot out the heatmap and or dotplot for GSEA results
# compile all the GSEA results from everyone 
cell_types = list.dirs(file.path("results", ANALYSIS_VERSION, 'wt13_enrichment'), recursive = FALSE, full.names = FALSE)
total_GSEA_results = data.frame()
for(cell_type in cell_types) {
  print(cell_type)
  GSEA_results = read.csv(file.path("results", ANALYSIS_VERSION, 'wt13_enrichment', cell_type, 'gsea_results_wt.csv'), row.names = 1)
  GSEA_results$cell_type = cell_type
  print(nrow(GSEA_results))
  if(nrow(total_GSEA_results) == 0) {
    total_GSEA_results = GSEA_results
  }
  else {
    total_GSEA_results = rbind(total_GSEA_results, GSEA_results)
  }
}

dir.create(file.path(TARGET_dir, "GSEA_heatmap"))

# create giant data 
interest_GSEA = c("cytoplasmic translation (GO:0002181)", 
                  "mRNA splicing, via spliceosome (GO:0000398)", 
                  "negative regulation of translation (GO:0017148)", 
                  "stem cell fate determination (GO:0048867)", 
                  "central nervous system development (GO:0007417)")
withr::with_dir(file.path(TARGET_dir, "GSEA_heatmap"), {
  interest_GSEA = vector()
  for(cell_type in unique(total_GSEA_results$cell_type)) {
    print(cell_type)
    
    temp_GSEA = total_GSEA_results[total_GSEA_results$cell_type == cell_type, ]
    temp_GSEA = temp_GSEA[temp_GSEA$padj < 0.05, ]
    temp_GSEA = temp_GSEA[temp_GSEA$NES > 0, ]
    temp_GSEA = temp_GSEA[order(temp_GSEA$NES, decreasing = TRUE), ]
    temp_GSEA = temp_GSEA[1:min(nrow(temp_GSEA), 3), ]
    temp_GSEA = temp_GSEA[order(temp_GSEA$NES, decreasing = FALSE), ]
    interest_GSEA = c(interest_GSEA, as.vector(temp_GSEA$pathway))
  }
  
  plot_df = total_GSEA_results[total_GSEA_results$pathway %in% interest_GSEA, ]
  plot_df[plot_df$padj > 0.05, 'NES'] = 0
  
  min_df = plot_df[, c('cell_type', 'pathway', 'NES')]
  
  df_wide = dcast(min_df, 
                  cell_type ~ pathway, 
                  value.var = "NES")
  rownames(df_wide) = df_wide$cell_type
  df_wide$cell_type = NULL
  cell_type_label = rownames(df_wide)
  df_wide = as.matrix(sapply(df_wide, as.numeric))  
  rownames(df_wide) = cell_type_label
  
  colclus <- hclust(dist(df_wide)) #cluster the columns
  new_order = colclus[["labels"]][colclus[['order']]]
  plot_df$cell_type = factor(plot_df$cell_type,levels=as.vector(new_order))
  
  p = ggplot(plot_df, aes(cell_type, pathway, fill= NES)) + 
    geom_tile() + 
    scale_x_discrete(position = "top") + 
    xlab("") +
    ggtitle(cell_type) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0)) + 
    scale_fill_gradient2(low="blue", mid = 'white', high="red", na.value="grey", midpoint = 0)
  
})

withr::with_dir(file.path(TARGET_dir, "GSEA_heatmap"), {
  for(cell_type in unique(total_GSEA_results$cell_type)) {
    print(cell_type)
    temp_GSEA = total_GSEA_results[total_GSEA_results$cell_type == cell_type, ]
    temp_GSEA = temp_GSEA[temp_GSEA$padj < 0.05, ]
    temp_GSEA = temp_GSEA[temp_GSEA$NES > 0, ]
    temp_GSEA = temp_GSEA[order(temp_GSEA$NES, decreasing = TRUE), ]
    temp_GSEA = temp_GSEA[1:min(nrow(temp_GSEA), 20), ]
    temp_GSEA = temp_GSEA[order(temp_GSEA$NES, decreasing = FALSE), ]
    interest_GSEA = as.vector(temp_GSEA$pathway)
    
    plot_df = total_GSEA_results[total_GSEA_results$pathway %in% interest_GSEA, ]
    plot_df$log_adj_p = -log(plot_df$padj)
    plot_df$pathway = factor(plot_df$pathway,levels=as.vector(temp_GSEA$pathway))
    plot_df[plot_df$padj > 0.05, 'NES'] = 0
    min_df = plot_df[, c('cell_type', 'pathway', 'NES')]
    
    df_wide = dcast(min_df, 
                    cell_type ~ pathway, 
                    value.var = "NES")
    rownames(df_wide) = df_wide$cell_type
    df_wide$cell_type = NULL
    cell_type_label = rownames(df_wide)
    df_wide = as.matrix(sapply(df_wide, as.numeric))  
    rownames(df_wide) = cell_type_label
    
    colclus <- hclust(dist(df_wide)) #cluster the columns
    new_order = colclus[["labels"]][colclus[['order']]]
    plot_df$cell_type = factor(plot_df$cell_type,levels=as.vector(new_order))
    
    p = ggplot(plot_df, aes(cell_type, pathway, fill= NES)) + 
      geom_tile() + 
      scale_x_discrete(position = "top") + 
      xlab("") +
      ggtitle(cell_type) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0)) + 
      scale_fill_gradient2(low="blue", mid = 'white', high="red", na.value="grey", midpoint = 0)
    
    ggsave(filename = paste0(cell_type, "_heatmap.png"), device = 'png', plot = p, width = 12, height = 4 + (length(unique(plot_df$pathway)) * 0.1))
    
    p = ggplot(plot_df, aes(x = cell_type, y = pathway)) + 
      geom_point(aes(size = log_adj_p, color = NES)) +
      scale_size_continuous(range = c(4,8)) +
      theme_bw(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0)) + 
      scale_x_discrete(position = "top") + 
      xlab("") +
      ggtitle(cell_type) +
      scale_colour_gradient2(high="red", mid = 'white', low = 'blue', midpoint = 0) +
      ylab('Normalized Enrichment Scores') 
    ggsave(filename = paste0(cell_type, "_dotplot.png"), device = 'png', plot = p, width = 12, height = 4 + (length(unique(plot_df$pathway)) * 0.1))
  }
})

dir.create(file.path(TARGET_dir, "GSEA_barplots"))

compile_bar_GSEA <- function(total_GSEA_results, cat) {
  sub_GSEA_results = total_GSEA_results[total_GSEA_results$pathway == cat, ]
  sub_GSEA_results = sub_GSEA_results[sub_GSEA_results$padj < 0.05, ]
  
  sub_GSEA_results$abs_NES = abs(sub_GSEA_results$NES)
  sub_GSEA_results$direction = NA
  sub_GSEA_results[sub_GSEA_results$NES > 0 , 'direction'] = 'Up'
  sub_GSEA_results[sub_GSEA_results$NES < 0, 'direction'] = 'Down'
  
  return(sub_GSEA_results)
}
withr::with_dir(file.path(TARGET_dir, "GSEA_barplots"), {
  cat_list = c("mRNA splicing, via spliceosome (GO:0000398)", 
               "translation (GO:0006412)", 
               "mitochondrial translation (GO:0032543)", 
               "apical junction assembly (GO:0043297)", 
               "septate junction assembly (GO:0019991)")
  for(cat in cat_list) {
    sub_GSEA_results = compile_bar_GSEA(total_GSEA_results, cat)
    p<-ggplot(data=sub_GSEA_results, aes(x=reorder(cell_type, NES), y=NES, fill=direction)) +
      geom_bar(stat="identity") + coord_flip() + theme_bw() + xlab("Cell Types") + 
      ggtitle(cat) + 
      scale_fill_manual(values=c("steel blue", "salmon"))
    ggsave(filename = paste0(cat, "_barplot.png"), plot = p, width = 8, height = 6)
  } 
})
# plot out the proportion of data 
proportion_df = data.frame("cell_types" = names(table(wt_object@meta.data$manual_celltypes)), 
                           "number_cells" = as.vector(table(wt_object@meta.data$manual_celltypes)))
proportion_df$cell_proportion = proportion_df$number_cells / sum(proportion_df$number_cells)

p<-ggplot(data=proportion_df, aes(x=reorder(cell_types, cell_proportion), y=cell_proportion, fill = cell_types)) +
  geom_bar(stat="identity") + theme_bw() + coord_flip() + 
  ylab("Total Cell Proportion") + 
  xlab("Cell Types")
ggsave(filename = file.path(TARGET_dir, "cell_proportion", paste0("cell_proportion.png")), plot = p, width = 14, height = 6)

##########################
# plot out the new marker genes 
wt_object@meta.data = wt_object@meta.data[, grepl("ModuleScore1", colnames(wt_object@meta.data)) == FALSE]
marker_genes = read.csv(file.path('results', ANALYSIS_VERSION, 'manual_annotation_wt13/manualMarkerGenes_2.csv'))
marker_genes_list = list()

for(unique_ct in unique(marker_genes$annotation)) {
  sub_marker_genes = marker_genes[marker_genes$annotation == unique_ct, ]
  marker_genes_list[[trimws(unique_ct)]] = intersect(as.vector(sub_marker_genes$marker_genes), rownames(wt_object@assays$RNA@data))
  wt_object = AddModuleScore(
    object = wt_object,
    features = list(marker_genes_list[[trimws(unique_ct)]]),
    name = paste0(stringr::str_replace_all(trimws(unique_ct), pattern = " ", replacement = "_"), '_ModuleScore')
  )
}

# plot out the feature plot 
module_score_names = colnames(wt_object@meta.data)
module_score_names = module_score_names[grep("ModuleScore", module_score_names)]
for(module_score_name in module_score_names) {
  geneset_name = stringr::str_remove_all(string = module_score_name, pattern = "_ModuleScore1")
  p = VlnPlot(wt_object, features = module_score_name, group.by = 'manual_celltypes', pt.size = 0) + 
    ggtitle(geneset_name) +  
    theme(plot.margin = margin(10, 10, 10, 20))
  ggsave(filename = file.path(TARGET_dir, "celltype_score", paste0(geneset_name, "_vlnplot.png")), plot = p, width = 16, height = 6)
  
  p = FeaturePlot(wt_object, features = module_score_name) + 
    ggtitle(geneset_name) 
  ggsave(filename = file.path(TARGET_dir, "celltype_score", paste0(geneset_name, "_UMAP.png")), plot = p, width = 8, height = 6)
}

for(temp_index in rownames(marker_genes)) { 
  temp_ct = marker_genes[temp_index, 'annotation']
  temp_gene = marker_genes[temp_index, 'marker_genes']
  if(temp_gene %in% rownames(wt_object) == FALSE) {
    next()
  }
  p = VlnPlot(wt_object, features = temp_gene, group.by = 'manual_celltypes', pt.size = 0) + 
    ggtitle(temp_gene) +  
    theme(plot.margin = margin(10, 10, 10, 20))
  ggsave(filename = file.path(TARGET_dir, "gene_scores", paste0(temp_ct, "_", temp_gene, "_vlnplot.png")), plot = p, width = 16, height = 6)
  
  p = FeaturePlot(wt_object, features = temp_gene) + 
    ggtitle(temp_gene) 
  ggsave(filename = file.path(TARGET_dir, "gene_scores", paste0(temp_ct, "_", temp_gene, "_UMAP.png")), plot = p, width = 8, height = 6)
}

################

