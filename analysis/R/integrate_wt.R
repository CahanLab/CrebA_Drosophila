library(Seurat)
library(harmony)
library(magrittr)
library(ggplot2)
library(stringr)

SAMPLETYPE = "wt"
TARGET_dir = file.path("results", ANALYSIS_VERSION, "wt_integrated")
dir.create(TARGET_dir)

rep1 = readRDS(file.path("results", ANALYSIS_VERSION, paste0(SAMPLETYPE, "_rep1"), "object.Rdata"))
rep1@meta.data$batch = 'rep_1'

rep2 = readRDS(file.path("results", ANALYSIS_VERSION, paste0(SAMPLETYPE, "_rep2"), "object.Rdata"))
rep2@meta.data$batch = 'rep_2'

wt_1_labels = read.csv("results/v18/Bianca_wt_cluster/wt1_clusters.csv")
colnames(wt_1_labels) = c('cluster_number', 'cluster_identification')

wt_2_labels = read.csv("results/v18/Bianca_wt_cluster/wt2_clusters.csv")
colnames(wt_2_labels) = c('cluster_number', 'cluster_identification')

assign_ct <- function(seurat_obj, annotate_st) { 
  seurat_obj@meta.data$Bianca_CT = NA
  for(cluster in unique(annotate_st$cluster_number)) { 
    seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters == str_trim(cluster), 'Bianca_CT'] = annotate_st[annotate_st$cluster_number == cluster, 'cluster_identification']
  }
  return(seurat_obj)
}

rep1 = assign_ct(rep1, wt_1_labels)
rep2 = assign_ct(rep2, wt_2_labels)

merged_seurat = merge(rep1, rep2)

raw_object = Seurat::CreateSeuratObject(merged_seurat@assays$RNA@counts, project = 'merged_samples')
raw_object@meta.data = merged_seurat@meta.data

saveRDS(raw_object, file = file.path(TARGET_dir, paste0(SAMPLETYPE, "_raw_merged.rds")))

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
object %<>% FindClusters()
pdf(file = file.path(TARGET_dir, "naive_integration.pdf"))
DimPlot(object, group.by = "batch")
dev.off()

# harmony integration 
object %<>% harmony::RunHarmony("batch")
object %<>% FindNeighbors(reduction = "harmony", dims = 1:num_pc)
object %<>% RunUMAP(reduction = "harmony", dim = 1:num_pc)
object %<>% FindClusters()
pdf(file = file.path(TARGET_dir, "harmony_integration.pdf"))
DimPlot(object, group.by = "batch")
dev.off()

saveRDS(object, file = file.path(TARGET_dir, 'object.rds'))

object = readRDS(file.path(TARGET_dir, 'object.rds'))
# get the smaller cluster annotation 
object@meta.data$smaller_tentativeCellType = NULL
for(i in rownames(metadata)[1:2]) { 
  i = as.numeric(i)
  individual_object = readRDS(file.path("results", ANALYSIS_VERSION, metadata$sample[i], 'object.Rdata'))
  individual_object_meta = individual_object@meta.data
  rownames(individual_object_meta) = paste0(rownames(individual_object_meta), "_", i)
  
  individual_clusterAnnotation = read.csv(file.path("results", ANALYSIS_VERSION, metadata$sample[i], 'tentativeCellTypes_BDGP.csv'), row.names = 1)
  individual_object_meta$tentativeCellType = NULL
  for(temp_cluster in unique(individual_clusterAnnotation$cluster)) { 
    individual_object_meta[individual_object_meta$seurat_clusters == temp_cluster, 'tentativeCellType'] = individual_clusterAnnotation[individual_clusterAnnotation$cluster == temp_cluster, 'annotation']
  }
  
  object@meta.data[rownames(individual_object_meta), 'smaller_tentativeCellType'] = individual_object_meta$tentativeCellType
}

DimPlot(object, group.by = 'Bianca_CT')
withr::with_dir(
  file.path(TARGET_dir), 
  { 
    
    curate_plot_df <- function(int_seurat, target_cell = 'all') { 
      UMAP_coord = int_seurat@reductions$umap@cell.embeddings
      UMAP_coord = as.data.frame(UMAP_coord)
      UMAP_coord$target_cells = NA
      UMAP_coord$batch = int_seurat$batch
      if(target_cell == 'all') { 
        UMAP_coord$target_cells = int_seurat$smaller_tentativeCellType
      } else{
        UMAP_coord[int_seurat$smaller_tentativeCellType == target_cell, 'target_cells'] = target_cell
        UMAP_coord[is.na(UMAP_coord$target_cells), 'target_cells'] = 'other'
      }
      return(UMAP_coord)
    }
    
    dir.create("overlayPlots_UMAPS")
    for(unique_celltype in unique(object@meta.data$smaller_tentativeCellType)) {
      plot_df = curate_plot_df(object, unique_celltype)
      plot_df[plot_df$target_cells != 'other', 'target_cells'] = paste0(plot_df[plot_df$target_cells != 'other', 'batch'], "_", plot_df[plot_df$target_cells != 'other', 'target_cells'])
      p = ggplot(plot_df) + 
        geom_point(data = subset(plot_df, target_cells == 'other'), aes(x=UMAP_1, y=UMAP_2, color = target_cells, alpha = 0.8))+
        geom_point(data = subset(plot_df, target_cells != 'other'), aes(x=UMAP_1, y=UMAP_2, color = target_cells, alpha = 0.8))+ 
        scale_color_manual(values=RColorBrewer::brewer.pal(12, 'Set3')[-c(2)]) + 
        ggtitle(unique_celltype)+
        theme_bw()+
        theme(legend.text = element_text(size = 12))
      
      ggsave(file.path('overlayPlots_UMAPS', paste0(unique_celltype, ".png")), plot = p, device = 'png', width = 12, height = 8)
    }  
  }
)

withr::with_dir(
  file.path(TARGET_dir), 
  { 
    
    curate_plot_df <- function(int_seurat, target_cell = 'all') { 
      UMAP_coord = int_seurat@reductions$umap@cell.embeddings
      UMAP_coord = as.data.frame(UMAP_coord)
      UMAP_coord$target_cells = NA
      UMAP_coord$batch = int_seurat$batch
      if(target_cell == 'all') { 
        UMAP_coord$target_cells = int_seurat$Bianca_CT
      } else{
        UMAP_coord[int_seurat$Bianca_CT == target_cell, 'target_cells'] = target_cell
        UMAP_coord[is.na(UMAP_coord$target_cells), 'target_cells'] = 'other'
      }
      return(UMAP_coord)
    }
    
    dir.create("overlayPlots_Bianca_UMAPS")
    for(unique_celltype in unique(object@meta.data$Bianca_CT)) {
      plot_df = curate_plot_df(object, unique_celltype)
      plot_df[plot_df$target_cells != 'other', 'target_cells'] = paste0(plot_df[plot_df$target_cells != 'other', 'batch'], "_", plot_df[plot_df$target_cells != 'other', 'target_cells'])
      p = ggplot(plot_df) + 
        geom_point(data = subset(plot_df, target_cells == 'other'), aes(x=UMAP_1, y=UMAP_2, color = target_cells, alpha = 0.8))+
        geom_point(data = subset(plot_df, target_cells != 'other'), aes(x=UMAP_1, y=UMAP_2, color = target_cells, alpha = 0.8))+ 
        scale_color_manual(values=RColorBrewer::brewer.pal(12, 'Set3')[-c(2)]) + 
        ggtitle(unique_celltype)+
        theme_bw()+
        theme(legend.text = element_text(size = 12))
      
      ggsave(file.path('overlayPlots_Bianca_UMAPS', paste0(unique_celltype, ".png")), plot = p, device = 'png', width = 12, height = 8)
    }  
  }
)



#############################
diff_genes_list = list()
overlaping_cts = intersect(unique(object@meta.data[object@meta.data$batch == 'rep_1', 'smaller_tentativeCellType']), unique(object@meta.data[object@meta.data$batch == 'rep_2', 'smaller_tentativeCellType']))
for(ct in overlaping_cts) { 
  sub_object = subset(x = object, cells = rownames(object@meta.data)[object@meta.data$smaller_tentativeCellType == ct])
  diff_genes = FindMarkers(sub_object, group.by = 'batch', ident.1 = 'rep_1', ident.2 = 'rep_2')
  diff_genes_list[[ct]] = diff_genes
}

up_reg_genes = rep(0, nrow(object@assays$RNA@meta.features))
down_reg_genes = rep(0, nrow(object@assays$RNA@meta.features))
names(up_reg_genes) = rownames(object@assays$RNA@meta.features)
names(down_reg_genes) = rownames(object@assays$RNA@meta.features)

for(ct in overlaping_cts) { 
  print(ct)
  diff_genes = diff_genes_list[[ct]]
  diff_genes = diff_genes[diff_genes$p_val_adj < 0.05, ]
  pos_diff_genes = diff_genes[diff_genes$avg_log2FC > 1, ]
  down_diff_genes = diff_genes[diff_genes$avg_log2FC < -1, ]
  
  up_reg_genes[rownames(pos_diff_genes)] = up_reg_genes[rownames(pos_diff_genes)] + 1
  down_reg_genes[rownames(down_diff_genes)] = down_reg_genes[rownames(down_diff_genes)] + 1
  
}
withr::with_dir(
  file.path(TARGET_dir), 
  { 
    down_reg_genes = down_reg_genes[down_reg_genes == 13]
    up_reg_genes = up_reg_genes[up_reg_genes == 13]
    
    saveRDS(down_reg_genes, file = 'down_reg_genes.rds')
    saveRDS(up_reg_genes, file = 'up_reg_genes.rds')
    
    write.csv(down_reg_genes, file = 'down_reg_genes.csv')
  }
)

# I ran the GO enrichment analysis of the down regulated genes
# it seems like the common theme is that muscle development is lacking in wt_rep1 



