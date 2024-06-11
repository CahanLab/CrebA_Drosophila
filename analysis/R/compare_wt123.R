library(Seurat)
library(harmony)
library(magrittr)
library(ggplot2)
library(stringr)

SAMPLETYPE = "wt"
TARGET_dir = file.path("results", ANALYSIS_VERSION, "wt_123_compare")
dir.create(TARGET_dir)

runs_list = list()
runs_list[[1]] = c('rep1', 'rep3')
runs_list[[2]] = c('rep2', 'rep3')
runs_list[[3]] = c('rep1', 'rep2')

curate_plot_df <- function(int_seurat, target_cell = 'all') { 
  UMAP_coord = int_seurat@reductions$umap@cell.embeddings
  UMAP_coord = as.data.frame(UMAP_coord)
  UMAP_coord$target_cells = NA
  UMAP_coord$batch = int_seurat$batch
  if(target_cell == 'all') { 
    UMAP_coord$target_cells = int_seurat$tentativeCellType
  } else{
    UMAP_coord[int_seurat$tentativeCellType == target_cell, 'target_cells'] = target_cell
    UMAP_coord[is.na(UMAP_coord$target_cells), 'target_cells'] = 'other'
  }
  return(UMAP_coord)
}

for(run in runs_list) {
  print(run)
  
  dir.create(file.path(TARGET_dir, paste0(run, collapse = '_')))
  wt_1 = readRDS(file.path("results", ANALYSIS_VERSION, paste0(SAMPLETYPE, "_", run[1]), "BDGP_automated_annotation_object.rds"))
  wt_1@meta.data$batch = run[1]
  
  wt_2 = readRDS(file.path("results", ANALYSIS_VERSION, paste0(SAMPLETYPE, "_", run[2]), "BDGP_automated_annotation_object.rds"))
  wt_2@meta.data$batch = run[2]
  
  merged_seurat = merge(wt_1, wt_2)
  
  raw_object = Seurat::CreateSeuratObject(merged_seurat@assays$RNA@counts, project = 'merged_samples')
  raw_object@meta.data = merged_seurat@meta.data
  
  object = Seurat::NormalizeData(raw_object)
  cellCycleMarkers = read.csv("accessory_data/cellCycleMarkers.csv", skip = 1, header = T)
  object %<>% CellCycleScoring(s.features = cellCycleMarkers$S.phase.markers., g2m.features = cellCycleMarkers$G2.M.phase.markers.)
  
  genes_in_object = rownames(object@assays$RNA@meta.features)
  
  mitochondrially_encoded_genes = read.table("accessory_data/mitochondrially_encoded_genes.tsv", header = F)[[1]]
  ribosomal_rna_genes           = read.table("accessory_data/ribosomal_rna_genes.tsv", header = F)[[1]]
  ribosomal_protein_genes       = grep("^rpl|^rps", ignore.case = T, value = T, genes_in_object)
  object[["mitochondrial_transcript_total_expression"]] =
    GetAssayData(object, "counts") %>%
    extract( convert_fbgn_to_symbol( mitochondrially_encoded_genes ), ) %>%
    colSums
  object[["ribosomal_rna_total_expression"]] =
    GetAssayData(object, "counts") %>%
    extract( convert_fbgn_to_symbol( ribosomal_rna_genes ), ) %>%
    colSums
  object[["ribosomal_protein_total_expression"]] =
    GetAssayData(object, "counts") %>%
    extract( ribosomal_protein_genes, ) %>%
    colSums
  object[["mitochondrial_transcript_norm_expression"]] = object[["mitochondrial_transcript_total_expression"]] / object$nCount_RNA
  object[["ribosomal_rna_norm_expression"           ]] = object[["ribosomal_rna_total_expression"           ]] / object$nCount_RNA
  object[["ribosomal_protein_norm_expression"       ]] = object[["ribosomal_protein_total_expression"       ]] / object$nCount_RNA
  object[["log10_nCount_RNA"]] = object[["nCount_RNA"]] %>% log10
  
  object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  object %<>% ScaleData(features = VariableFeatures(object = object))
  object %<>% RunPCA(features = VariableFeatures(object = object), npcs = 100)
  pc_cutoff = RMTstat::qmp( 1, ndf=length(Cells(object)), pdim=length(VariableFeatures(object)), var=1)
  singular_values = slot(Reductions(object, slot = "pca"), "stdev")
  is_significant = singular_values^2 > pc_cutoff
  num_pc = sum(is_significant)
  
  # after getting the significant PCs, rerun PCA 
  object %<>% RunPCA(features = VariableFeatures(object = object), npcs = num_pc)
  object %<>% FindNeighbors(dims = 1:num_pc)
  object %<>% RunUMAP(dim = 1:num_pc)
  #object %<>% FindClusters()
  
  
  pdf(file = file.path(TARGET_dir, paste0(run, collapse = '_'), "naive_integration.pdf"))
  DimPlot(object, group.by = "batch")
  dev.off()
  
  # harmony integration 
  object %<>% harmony::RunHarmony("batch")
  object %<>% FindNeighbors(reduction = "harmony", dims = 1:num_pc)
  object %<>% RunUMAP(reduction = "harmony", dim = 1:num_pc)
  object %<>% FindClusters()
  pdf(file = file.path(TARGET_dir, paste0(run, collapse = '_'), "harmony_integration.pdf"))
  DimPlot(object, group.by = "batch")
  dev.off()
  
  for(ct in unique(object$tentativeCellType)) {
    plot_df = curate_plot_df(int_seurat = object, target_cell = ct)
    plot_df[plot_df$target_cells != 'other', 'target_cells'] = paste0(plot_df[plot_df$target_cells != 'other', 'batch'], "_", plot_df[plot_df$target_cells != 'other', 'target_cells'])
    
    p = ggplot(plot_df) + 
      geom_point(data = subset(plot_df, target_cells == 'other'), aes(x=UMAP_1, y=UMAP_2, color = target_cells, alpha = 0.8))+
      geom_point(data = subset(plot_df, target_cells != 'other'), aes(x=UMAP_1, y=UMAP_2, color = target_cells, alpha = 0.8))+ 
      scale_color_manual(values=RColorBrewer::brewer.pal(12, 'Set3')[-c(2)]) + 
      ggtitle(ct)+
      theme_bw()+
      theme(legend.text = element_text(size = 12))
    ggsave(file.path(TARGET_dir, paste0(run, collapse = '_'), paste0(ct, "_overlay.png")), plot = p, width = 8, height = 6)
    print(ct)
  }
  saveRDS(object, file = file.path(TARGET_dir, paste0(run, collapse = '_'), 'harmony_object.rds'))
}

# get the proportion plot dataframe 
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
# add it back on 
for(run in runs_list) {
  object = readRDS(file = file.path(TARGET_dir, paste0(run, collapse = '_'), 'harmony_object.rds'))
  proportion_plot = proportion_study(object)
  proportion_plot = proportion_plot[!duplicated(proportion_plot), ]
  p = ggplot(proportion_plot, aes(x=cluster_id, y=proportion, fill=batch)) +
    geom_bar(stat="identity", alpha=0.7) +
    theme_bw()
  ggsave(file.path(TARGET_dir, paste0(run, collapse = '_'), 'cluster_proportion.png'), plot = p, width = 10, height = 6)
  FeaturePlot(object, "log10_nCount_RNA")
  FeaturePlot(object, "nFeature_RNA")
  DimPlot(object, group.by = "batch")
  
  VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "ribosomal_protein_total_expression"), pt.size = 0, ncol = 3, group.by = 'batch')
}

##############
# merge 3 and see how that works 
wt_1 = readRDS(file.path("results", ANALYSIS_VERSION, paste0(SAMPLETYPE, "_rep1"), "object.Rdata"))
wt_1@meta.data$batch = "rep1"

wt_2 = readRDS(file.path("results", ANALYSIS_VERSION, paste0(SAMPLETYPE, "_rep2"), "object.Rdata"))
wt_2@meta.data$batch = 'rep2'

wt_3 = readRDS(file.path("results", ANALYSIS_VERSION, paste0(SAMPLETYPE, "_rep3"), "object.Rdata"))
wt_3@meta.data$batch = 'rep3'

merged_seurat = merge(wt_1, c(wt_2, wt_3), )
dir.create(file.path(TARGET_dir, paste0("all_3")))

raw_object = Seurat::CreateSeuratObject(merged_seurat@assays$RNA@counts, project = 'merged_samples')
raw_object@meta.data = merged_seurat@meta.data

object = Seurat::NormalizeData(raw_object)
cellCycleMarkers = read.csv("accessory_data/cellCycleMarkers.csv", skip = 1, header = T)
object %<>% CellCycleScoring(s.features = cellCycleMarkers$S.phase.markers., g2m.features = cellCycleMarkers$G2.M.phase.markers.)

genes_in_object = rownames(object@assays$RNA@meta.features)

mitochondrially_encoded_genes = read.table("accessory_data/mitochondrially_encoded_genes.tsv", header = F)[[1]]
ribosomal_rna_genes           = read.table("accessory_data/ribosomal_rna_genes.tsv", header = F)[[1]]
ribosomal_protein_genes       = grep("^rpl|^rps", ignore.case = T, value = T, genes_in_object)
object[["mitochondrial_transcript_total_expression"]] =
  GetAssayData(object, "counts") %>%
  extract( convert_fbgn_to_symbol( mitochondrially_encoded_genes ), ) %>%
  colSums
object[["ribosomal_rna_total_expression"]] =
  GetAssayData(object, "counts") %>%
  extract( convert_fbgn_to_symbol( ribosomal_rna_genes ), ) %>%
  colSums
object[["ribosomal_protein_total_expression"]] =
  GetAssayData(object, "counts") %>%
  extract( ribosomal_protein_genes, ) %>%
  colSums
object[["mitochondrial_transcript_norm_expression"]] = object[["mitochondrial_transcript_total_expression"]] / object$nCount_RNA
object[["ribosomal_rna_norm_expression"           ]] = object[["ribosomal_rna_total_expression"           ]] / object$nCount_RNA
object[["ribosomal_protein_norm_expression"       ]] = object[["ribosomal_protein_total_expression"       ]] / object$nCount_RNA
object[["log10_nCount_RNA"]] = object[["nCount_RNA"]] %>% log10

object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
object %<>% ScaleData(features = VariableFeatures(object = object))
object %<>% RunPCA(features = VariableFeatures(object = object), npcs = 100)
pc_cutoff = RMTstat::qmp( 1, ndf=length(Cells(object)), pdim=length(VariableFeatures(object)), var=1)
singular_values = slot(Reductions(object, slot = "pca"), "stdev")
is_significant = singular_values^2 > pc_cutoff
num_pc = sum(is_significant)

# after getting the significant PCs, rerun PCA 
object %<>% RunPCA(features = VariableFeatures(object = object), npcs = num_pc)
object %<>% FindNeighbors(dims = 1:num_pc)
object %<>% RunUMAP(dim = 1:num_pc)
#object %<>% FindClusters()


pdf(file = file.path(TARGET_dir, paste0("all_3"), "naive_integration.pdf"))
DimPlot(object, group.by = "batch")
dev.off()

# harmony integration 
object %<>% harmony::RunHarmony("batch")
object %<>% FindNeighbors(reduction = "harmony", dims = 1:num_pc)
object %<>% RunUMAP(reduction = "harmony", dim = 1:num_pc)
object %<>% FindClusters()
pdf(file = file.path(TARGET_dir, paste0("all_3"), "harmony_integration.pdf"))
DimPlot(object, group.by = "batch")
dev.off()

saveRDS(object, file = file.path(TARGET_dir, paste0("all_3"), 'harmony_object.rds'))

#################
library(ggplot2)
object = readRDS(file.path(TARGET_dir, paste0("all_3"), 'harmony_object.rds'))
VlnPlot(object, features = 'log10_nCount_RNA', group.by = 'batch', pt.size = 0)
ggsave(file.path(TARGET_dir, paste0("all_3"), 'log10_nCounts_RNA.png'), width = 8, height = 6)

VlnPlot(object, features = 'nFeature_RNA', group.by = 'batch', pt.size = 0)
ggsave(file.path(TARGET_dir, paste0("all_3"), 'nFeature_RNA.png'), width = 8, height = 6)

VlnPlot(object, features = 'mitochondrial_transcript_total_expression', group.by = 'batch', pt.size = 0)
ggsave(file.path(TARGET_dir, paste0("all_3"), 'mitochondrial_transcript_total_expression.png'), width = 8, height = 6)
