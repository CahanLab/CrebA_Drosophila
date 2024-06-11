library(Seurat)
library(stringr)
library(Matrix)
library(pheatmap)
library(enrichR)
library(RColorBrewer)
library(openxlsx)

cell_type = 'tracheal_system_wt123'
TARGET_dir = file.path("results", ANALYSIS_VERSION, cell_type)

wt_1 = readRDS(file.path("results", ANALYSIS_VERSION, paste0('wt', "_rep1"), "object.Rdata"))
wt_1@meta.data$batch = "rep1"
wt_1_ct = read.csv(file.path("results", ANALYSIS_VERSION, paste0('wt', "_rep1"), "tentativeCellTypes_BDGP.csv"), row.names = 1)
cluster_id = wt_1_ct[wt_1_ct$annotation == 'embryonic/larval tracheal system', 'cluster']
wt_1 = subset(wt_1, subset = seurat_clusters == cluster_id)

wt_2 = readRDS(file.path("results", ANALYSIS_VERSION, paste0('wt', "_rep2"), "object.Rdata"))
wt_2@meta.data$batch = 'rep2'
wt_2_ct = read.csv(file.path("results", ANALYSIS_VERSION, paste0('wt', "_rep2"), "tentativeCellTypes_BDGP.csv"), row.names = 1)
cluster_id = wt_2_ct[wt_2_ct$annotation == 'embryonic/larval tracheal system', 'cluster']
wt_2 = subset(wt_2, subset = seurat_clusters == cluster_id)

wt_3 = readRDS(file.path("results", ANALYSIS_VERSION, paste0('wt', "_rep3"), "object.Rdata"))
wt_3@meta.data$batch = 'rep3'
wt_3_ct = read.csv(file.path("results", ANALYSIS_VERSION, paste0('wt', "_rep3"), "tentativeCellTypes_BDGP.csv"), row.names = 1)
cluster_id = wt_3_ct[wt_3_ct$annotation == 'embryonic/larval tracheal system', 'cluster']
wt_3 = subset(wt_3, subset = seurat_clusters == cluster_id)

merged_seurat = merge(wt_1, c(wt_2, wt_3))

raw_object = Seurat::CreateSeuratObject(merged_seurat@assays$RNA@counts, project = 'merged_samples')
raw_object@meta.data = merged_seurat@meta.data

object = Seurat::NormalizeData(raw_object)

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
ggsave(filename = file.path(TARGET_dir, 'naive_cluster.png'), width = 8, height = 6)

DimPlot(object, group.by = 'batch')
ggsave(filename = file.path(TARGET_dir, 'naive_experimental_condition.png'), width = 8, height = 6)

DimPlot(object, group.by = 'batch')
ggsave(filename = file.path(TARGET_dir, 'naive_batch.png'), width = 8, height = 6)

saveRDS(object, file = file.path(TARGET_dir, "ct_object.rds"))
############
object %<>% harmony::RunHarmony("batch")
object %<>% FindNeighbors(reduction = "harmony", dims = 1:num_pc)
object %<>% RunUMAP(reduction = "harmony", dim = 1:num_pc)
object %<>% FindClusters()

DimPlot(object)
ggsave(filename = file.path(TARGET_dir, 'harmony_cluster.png'), width = 8, height = 6)

DimPlot(object, group.by = 'batch')
ggsave(filename = file.path(TARGET_dir, 'harmony_batch.png'), width = 8, height = 6)

cell_type = 'tracheal_system'
dir.create(file.path(TARGET_dir, "time_genes"))
TARGET_dir = file.path(TARGET_dir, "time_genes")
excel_path = file.path("accessory_data/continuum_drosophila_embryonic_development/Andrew_Marker_Genes", paste0(cell_type, "_time_genes.xlsx"))
page_names = openxlsx::getSheetNames(excel_path)

marker_genes_list = list()
union_genes = c()
all_genes = c()
for(page_name in page_names) {
  print(page_name)
  time_name = stringr::str_split(page_name, " hr")[[1]][1]
  time_name = paste0(" hr", time_name)
  time_name = stringr::str_remove_all(time_name, " ")
  time_name = stringr::str_replace_all(time_name, "-", "_")
  marker_page = openxlsx::read.xlsx(excel_path, sheet = page_name, colNames = FALSE)
  colnames(marker_page) = c("gene", 
                            "p_val", 
                            "avg_log2FC", 
                            "pct.1", 
                            "pct.2", 
                            "p_val_adj", 
                            "cluster", 
                            "time")
  rownames(marker_page) = marker_page$gene
  marker_page = marker_page[intersect(marker_page$gene, rownames(object@assays$RNA@meta.features)), ]
  marker_genes_list[[time_name]] = marker_page
  
  if(length(union_genes) == 0) { 
    union_genes = marker_page$gene
  }
  else {
    union_genes = intersect(union_genes, marker_page$gene)
  }
  
  all_genes = c(all_genes, marker_page$gene)
  object <- AddModuleScore(
    object = object,
    features = list(marker_page$gene),
    name = time_name
  )
  
  colnames(object@meta.data) <-
    gsub(x = colnames(object@meta.data)
         , pattern = paste0(time_name,1)
         , replacement = time_name
    )
  
  p = FeaturePlot(object, features = time_name)
  ggsave(filename = file.path(TARGET_dir, paste0(time_name, "_UMAP.png")), plot = p, width = 8, height = 6)
  
  p = VlnPlot(object, features = time_name, group.by = 'batch', pt.size = 0)
  ggsave(filename = file.path(TARGET_dir, paste0(time_name, "_vln.png")), plot = p, width = 8, height = 6)
}

# plotout the timepoint specific genes 
unique_genes = names(table(all_genes)[table(all_genes) == 1])

for(time_frame in names(marker_genes_list)) {
  print(time_frame)
  
  marker_page = marker_genes_list[[time_frame]]
  marker_page = marker_page[marker_page$gene %in% unique_genes, ]
  
  if(nrow(marker_page) == 0) {
    next
  }
  else{
    object <- AddModuleScore(
      object = object,
      features = list(marker_page$gene),
      name = time_name
    )
    
    colnames(object@meta.data) <-
      gsub(x = colnames(object@meta.data)
           , pattern = paste0(time_name,1)
           , replacement = time_name
      )
  }
  
  p = FeaturePlot(object, features = time_name)
  ggsave(filename = file.path(TARGET_dir, paste0(time_name, "_specific_UMAP.png")), plot = p, width = 8, height = 6)
  
  p = VlnPlot(object, features = time_name, group.by = 'batch', pt.size = 0)
  ggsave(filename = file.path(TARGET_dir, paste0(time_name, "_specific_vln.png")), plot = p, width = 8, height = 6)
}

#object = subset(object, subset = batch != 'rep3')
DE_markers = Seurat::FindMarkers(object = object, ident.1 = 'rep1', group.by = 'batch', test.use = 'bimod')
DE_markers = DE_markers[DE_markers$p_val_adj < 0.05, ]
rep2_genes = DE_markers[DE_markers$avg_log2FC > 0, ]

enrichment_results = enrichR::enrichr(
  genes = rownames(rep2_genes), 
  databases = c(
    "GO_Biological_Process_2018"
  )
)
biological_analysis = enrichment_results$GO_Biological_Process_2018
biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
