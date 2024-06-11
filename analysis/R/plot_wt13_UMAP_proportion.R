library(ggplot2)
library(RColorBrewer)
library(ggdendroplot)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "figure_plots", 'wt_late_UMAP_proportion')
dir.create(TARGET_dir, recursive = TRUE)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
qual_col_pals = qual_col_pals[rownames(qual_col_pals) != 'Dark2', ]
qual_col_pals = qual_col_pals[rownames(qual_col_pals) != 'Set2', ]

col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = col_vector[col_vector != col_vector[4]]
col_vector = col_vector[col_vector != col_vector[4]]

wt_object = readRDS(file.path('results', ANALYSIS_VERSION, 'manual_annotation_wt13/manual_celltype_object3.rds'))
wt_object@meta.data[wt_object@meta.data$batch == 'rep_3', 'batch'] = 'rep_2'

withr::with_dir(TARGET_dir, {
  p = DimPlot(wt_object, group.by = 'manual_celltypes', label = TRUE, label.size = 5) +
    ggtitle("Celltype Labels") + 
    xlim(c(-8, 12))
  ggsave(file.path("cell_type_UMAP.png"), plot = p, width = 20, height = 10)
  
  p = DimPlot(wt_object, group.by = 'batch', label = FALSE, label.size = 5) +
    ggtitle("Batch") + 
    xlim(c(-8, 12))
  ggsave(file.path("Batch_UMAP.png"), plot = p, width = 20, height = 10)
  
  p = DimPlot(wt_object, group.by = 'seurat_clusters', label = TRUE, label.size = 8) +
    ggtitle("Seurat Clusters") + 
    xlim(c(-8, 12))
  ggsave(file.path("clusters_UMAP.png"), plot = p, width = 20, height = 10)
  
})

# plot out the proportion of each cell type 
withr::with_dir(TARGET_dir, {
  proportion_df = data.frame("cell_types" = names(table(wt_object@meta.data$manual_celltypes)), 
                             "number_cells" = as.vector(table(wt_object@meta.data$manual_celltypes)))
  proportion_df$cell_proportion = proportion_df$number_cells / sum(proportion_df$number_cells)
  
  p<-ggplot(data=proportion_df, aes(x=reorder(cell_types, cell_proportion), y=cell_proportion, fill = cell_types)) +
    scale_fill_discrete(name = "Cell Types") +
    geom_bar(stat="identity") + theme_bw() + coord_flip() + 
    ylab("Total Cell Proportion") + 
    xlab("Cell Types")
  ggsave(filename = file.path("cell_proportion_bar.png"), plot = p, width = 14, height = 6)
})

# plot out the proportion of each cell type per batch
withr::with_dir(TARGET_dir, {
  proportion_df = data.frame()
  for(batch in unique(wt_object@meta.data$batch)) {
    sub_wt_object = subset(wt_object, subset = batch == batch)
    temp_proportion_df = data.frame("cell_types" = names(table(sub_wt_object@meta.data$manual_celltypes)), 
                               "number_cells" = as.vector(table(sub_wt_object@meta.data$manual_celltypes)))
    temp_proportion_df$cell_proportion = temp_proportion_df$number_cells / sum(temp_proportion_df$number_cells)
    temp_proportion_df$batch = batch
    proportion_df = rbind(proportion_df, temp_proportion_df)
  }

  
  p<-ggplot(data=proportion_df, aes(x=reorder(cell_types, cell_proportion), y=cell_proportion, fill = batch)) +
    scale_fill_discrete(name = "Cell Types") +
    geom_bar(stat="identity", position = 'dodge') + theme_bw() + coord_flip() + 
    ylab("Total Cell Proportion") + 
    xlab("Cell Types")
  ggsave(filename = file.path("cell_proportion_bar.png"), plot = p, width = 14, height = 6)
})
