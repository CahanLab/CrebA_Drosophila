library(stringr)
library(ggplot2)
library(bluster)
library(cluster)
library(RColorBrewer)
library(plyr)
#TODO -- write the rough draft script to fit this 
#TODO -- transplant my old script to this 
##########################################

##########################################
proportion_study <- function(object) { 
  plot_df = matrix(nrow = 2, ncol = 3)
  colnames(plot_df) = c('cluster_id', 'batch', 'proportion')
  plot_df = as.data.frame(plot_df)
  plot_df$cluster_id = 'full_data'
  plot_df$batch = names(table(object$batch))
  plot_df$proportion = table(object$batch) / length(object$batch)
  
  meta_tab = object@meta.data
  for(temp_cluster in unique(object$seurat_clusters)) { 
    temp_meta = meta_tab[meta_tab$seurat_clusters == temp_cluster, ]
    temp_plot_df = matrix(nrow = 2, ncol = 3)
    colnames(temp_plot_df) = c('cluster_id', 'batch', 'proportion')
    temp_plot_df = as.data.frame(temp_plot_df)
    temp_plot_df$cluster_id = temp_cluster
    temp_plot_df$batch = names(table(temp_meta$batch))
    temp_plot_df$proportion = table(temp_meta$batch) / length(temp_meta$batch)
    
    plot_df = rbind(plot_df, temp_plot_df)
  }
  return(plot_df)
}

integration_methods = list.dirs(file.path("results", ANALYSIS_VERSION, "test_batch_correction_100pc"), full.names = FALSE, recursive = FALSE)

mean_diff = c()
std_diff = c()
for(integration_method in integration_methods) {
  object = readRDS(file.path("results", ANALYSIS_VERSION, "test_batch_correction_100pc", integration_method, 'processed_object.rds'))
  if(integration_method == 'liger') { 
    str_split_df = stringr::str_split_fixed(rownames(object@meta.data), "_", n = 3)
    object@meta.data$batch = paste0(str_split_df[, 1], "_", str_split_df[, 2])
  }
  prop_plot_df = proportion_study(object)
  prop_plot_df = prop_plot_df[prop_plot_df$batch == 'rib_1', ]
  rownames(prop_plot_df) = prop_plot_df$cluster_id
  
  prop_plot_df$proportion = abs(prop_plot_df$proportion - prop_plot_df['full_data', 'proportion'])
  prop_plot_df = prop_plot_df[rownames(prop_plot_df) != 'full_data', ]
  
  mean_diff = c(mean_diff, mean(prop_plot_df$proportion))
  std_diff = c(std_diff, sd(prop_plot_df$proportion))
}

plot_df = data.frame('integration_method' = integration_methods, 
                     'mean_diff' = mean_diff, 
                     'stdv_diff' = std_diff)

p = ggplot(plot_df) +
  geom_bar( aes(x=integration_method, y=mean_diff), stat="identity", fill="skyblue", alpha=0.7) +
  geom_errorbar( aes(x=integration_method, ymin=mean_diff-stdv_diff, ymax=mean_diff+stdv_diff), width=0.4, colour="orange", alpha=0.9, size=1.3) + 
  theme_bw()
ggsave(file.path("results", ANALYSIS_VERSION, "test_batch_correction_100pc", 'mixing_level.png'), plot = p, )
##########################################
seurat_rib1 = readRDS(file.path("results", ANALYSIS_VERSION, "rib_rep1/object.Rdata"))
meta_rib1 = seurat_rib1@meta.data
meta_rib1$cell_type = NA
ct_rib1 = read.csv(file.path("results", ANALYSIS_VERSION, 'rib_rep1/tentativeCellTypes.csv'), row.names = 1)
for(cell_cluster in ct_rib1$cluster) { 
  meta_rib1[meta_rib1$seurat_clusters == cell_cluster, 'cell_type'] = ct_rib1[ct_rib1$cluster == cell_cluster, 'annotation']
}
seurat_rib1@meta.data = meta_rib1
seurat_rib1@meta.data$batch = 'rib_1'

seurat_rib2 = readRDS(file.path("results", ANALYSIS_VERSION, "rib_rep2/object.Rdata"))
meta_rib2 = seurat_rib2@meta.data
meta_rib2$cell_type = NA
ct_rib2 = read.csv(file.path("results", ANALYSIS_VERSION, 'rib_rep2/tentativeCellTypes.csv'), row.names = 1)
for(cell_cluster in ct_rib2$cluster) { 
  meta_rib2[meta_rib2$seurat_clusters == cell_cluster, 'cell_type'] = ct_rib2[ct_rib2$cluster == cell_cluster, 'annotation']
}
seurat_rib2@meta.data = meta_rib2
seurat_rib2@meta.data$batch = 'rib_2'

# check for intersecting cells 
i_cells = intersect(unique(seurat_rib1$cell_type), unique(seurat_rib2$cell_type))

seurat_combined = merge(seurat_rib1, seurat_rib2)
meta_combined = seurat_combined@meta.data
meta_combined = meta_combined[!is.na(meta_combined$cell_type), ]
meta_combined = meta_combined[meta_combined$cell_type %in% i_cells, ]

integration_methods = list.dirs(file.path("results", ANALYSIS_VERSION, "test_batch_correction_100pc"), full.names = FALSE, recursive = FALSE)

curate_plot_df <- function(int_seurat, target_cell = 'all') { 
  UMAP_coord = int_seurat@reductions$umap@cell.embeddings
  UMAP_coord = as.data.frame(UMAP_coord)
  UMAP_coord$target_cells = NA
  
  if(target_cell == 'all') { 
    UMAP_coord$target_cells = int_seurat$cell_type
  } else{
    UMAP_coord[int_seurat$cell_type == target_cell & int_seurat$batch == 'rib_1', 'target_cells'] = paste0('rib_1_', target_cell)
    UMAP_coord[int_seurat$cell_type == target_cell & int_seurat$batch == 'rib_2', 'target_cells'] = paste0('rib_2_', target_cell)
    UMAP_coord[is.na(UMAP_coord$target_cells), 'target_cells'] = 'rib_other'
  }
  return(UMAP_coord)
}

for(integration_method in integration_methods) { 
  int_seurat = readRDS(file.path("results", ANALYSIS_VERSION, "test_batch_correction_100pc", integration_method, "processed_object.rds"))
  if(integration_method == 'liger') { # liger is special in the integration 
    int_PCA = int_seurat@reductions$inmf@cell.embeddings
    rownames(int_PCA) = stringr::str_split(rownames(int_PCA), pattern = "_", n = 3, simplify = TRUE)[, 3]
    int_seurat@meta.data$batch = paste0(stringr::str_split(rownames(int_seurat@meta.data), pattern = "_", n = 3, simplify = TRUE)[, 1], "_", stringr::str_split(rownames(int_seurat@meta.data), pattern = "_", n = 3, simplify = TRUE)[, 2])
    int_seurat = Seurat::RenameCells(int_seurat, new.names = stringr::str_split(rownames(int_seurat@meta.data), pattern = "_", n = 3, simplify = TRUE)[, 3])
  }
  else if(integration_method == 'naive'){ 
    int_PCA = int_seurat@reductions$pca@cell.embeddings
  }
  else if(integration_method == 'harmony'){ 
    int_PCA = int_seurat@reductions$harmony@cell.embeddings
  }
  else if(integration_method == 'seurat'){ 
    int_PCA = int_seurat@reductions$pca@cell.embeddings
  }
  
  int_PCA = int_PCA[rownames(meta_combined), ]
  int_PCA = as.data.frame(int_PCA)
  #int_PCA = int_PCA[, 1:20] # only use first 20 PCs as described in the benchmarking paper
  
  for_python_PCA = int_PCA
  for_python_PCA$cell_type = meta_combined$cell_type
  write.csv(for_python_PCA, file = file.path("results", ANALYSIS_VERSION, "test_batch_correction_100pc", integration_method, "for_python_pca.csv"))
  
  fast_Silhouette = bluster::approxSilhouette(int_PCA, meta_combined$cell_type)
  saveRDS(fast_Silhouette, file = file.path("results", ANALYSIS_VERSION, "test_batch_correction_100pc", integration_method, "cell_type_silhouette_analysis.rds"))
  int_seurat = subset(int_seurat, cells = colnames(seurat_combined))
  int_seurat@meta.data$cell_type = seurat_combined$cell_type

  plot_df = curate_plot_df(int_seurat)
  p = ggplot(plot_df, aes(x=UMAP_1, y=UMAP_2, color = target_cells)) + 
    geom_point() + 
    #scale_color_manual(values=RColorBrewer::brewer.pal(12, 'Paired')) + 
    theme_bw()
  
  dir.create(file.path("results", ANALYSIS_VERSION, "test_batch_correction_100pc", integration_method, "UMAP_celltypes/"))
  ggsave(file.path("results", ANALYSIS_VERSION, "test_batch_correction_100pc", integration_method, "UMAP_celltypes/all.png"), plot = p, device = 'png', width = 12, height = 8)
  
  for(unique_celltype in unique(meta_combined$cell_type)) {
    plot_df = curate_plot_df(int_seurat, unique_celltype)
    
    p = ggplot(plot_df) + 
      geom_point(data = subset(plot_df, target_cells == sort(unique(plot_df$target_cells))[3]), aes(x=UMAP_1, y=UMAP_2, color = target_cells))+
      geom_point(data = subset(plot_df, target_cells == sort(unique(plot_df$target_cells))[2]), aes(x=UMAP_1, y=UMAP_2, color = target_cells)) + 
      geom_point(data = subset(plot_df, target_cells == sort(unique(plot_df$target_cells))[1]), aes(x=UMAP_1, y=UMAP_2, color = target_cells)) + 
      scale_color_manual(values=RColorBrewer::brewer.pal(12, 'Set3')[c(4, 5, 9)]) + 
      theme_bw()
    ggsave(file.path("results", ANALYSIS_VERSION, "test_batch_correction_100pc", integration_method, "UMAP_celltypes", paste0(unique_celltype, ".png")), plot = p, device = 'png', width = 12, height = 8)
  }
}

mean_sil_vector = vector()
for(integration_method in integration_methods) { 
  fast_Silhouette = readRDS(file.path("results/v16/test_batch_correction_100pc", integration_method, "cell_type_silhouette_analysis.rds"))
  
  ct_mean_sil = vector()
  
  for(unique_ct in unique(fast_Silhouette$cluster)) {
    temp_Sil = fast_Silhouette[fast_Silhouette$cluster == unique_ct, ]
    ct_mean_sil = c(ct_mean_sil, mean(temp_Sil$width))
  }
  ct_sil_df = data.frame(unique(fast_Silhouette$cluster), ct_mean_sil)
  colnames(ct_sil_df) = c("cell_types", 'mean_silhouette')
  write.csv(ct_sil_df, file = file.path("results/v16/test_batch_correction_100pc", integration_method, "ct_sil_df.csv"))
  mean_sil_vector = c(mean_sil_vector, mean(ct_mean_sil))
  
}
mean_sil_df = data.frame(integration_methods, mean_sil_vector)
write.csv(mean_sil_df, file = file.path("results/v16/test_batch_correction_100pc", "mean_sil_df.csv"))

for(integration_method in integration_methods) { 
  ct_sil_df = read.csv(file.path("results/v16/test_batch_correction_100pc", integration_method, "ct_sil_df.csv"), row.names = 1)
  p = ggplot(ct_sil_df) +
    geom_bar( aes(x=cell_types, y=mean_silhouette), stat="identity", fill="skyblue", alpha=0.7) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(file.path("results", ANALYSIS_VERSION, "test_batch_correction_100pc", integration_method, 'average_silhouette.png'), plot = p, )
}

mean_sil_df = read.csv(file.path("results/v16/test_batch_correction_100pc", "mean_sil_df.csv"), row.names = 1)
p = ggplot(mean_sil_df) +
  geom_bar( aes(x=integration_methods, y=mean_sil_vector), stat="identity", fill="skyblue", alpha=0.7) +
  ylab("Average Silhouette") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave(file.path("results", ANALYSIS_VERSION, "test_batch_correction_100pc", 'average_silhouette_acrossCT.png'), plot = p, width = 7, height = 4)

