library(Seurat)
library(stringr)
library(singleCellNet)
library(Matrix)
library(pheatmap)
library(enrichR)
library(RColorBrewer)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "tracheal_system")
dir.create(TARGET_dir)

wt_object = readRDS(file.path("results", ANALYSIS_VERSION, "wt13_integrated", "BDGP_automated_annotation_object.rds"))
rib_object = readRDS(file.path("results", ANALYSIS_VERSION, "rib_integrated/", "BDGP_automated_annotation_object.rds"))

ts_wt_obj = subset(wt_object, subset = Integrated_tentativeCellType == 'tracheal system') # subset out the salivary gland
ts_wt_obj@meta.data$experimental_condition = 'Wt'

ts_rib_obj = subset(rib_object, subset = Integrated_tentativeCellType == 'tracheal system') # subset out the salivary gland
ts_rib_obj@meta.data$experimental_condition = 'Rib'

combined_obj = merge(ts_wt_obj, ts_rib_obj)
object = Seurat::CreateSeuratObject(combined_obj@assays$RNA@counts, project = 'ts_cells')
object@meta.data$experimental_condition = combined_obj@meta.data$experimental_condition
object@meta.data$batch = paste0(combined_obj@meta.data$experimental_condition, "_", combined_obj@meta.data$batch)

object = Seurat::NormalizeData(object)

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

DimPlot(object, group.by = 'experimental_condition')
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

DimPlot(object, group.by = 'experimental_condition')
ggsave(filename = file.path(TARGET_dir, 'harmony_experimental_condition.png'), width = 8, height = 6)

DimPlot(object, group.by = 'batch')
ggsave(filename = file.path(TARGET_dir, 'harmony_batch.png'), width = 8, height = 6)

object = readRDS(file.path(TARGET_dir, "ct_object.rds"))
markers <- FindMarkers(object = object, ident.1 = "Wt", group.by = 'experimental_condition', test.use = 'wilcox') 
#markers = markers[markers$avg_log2FC < 0, ]
write.csv(markers, file = file.path(TARGET_dir, 'DE_genes_wt_rib_ts_wilcox.csv'))

markers <- FindMarkers(object = object, ident.1 = "Wt", group.by = 'experimental_condition', test.use = 'bimod') 

#markers = markers[markers$avg_log2FC < 0, ]
write.csv(markers, file = file.path(TARGET_dir, 'DE_genes_wt_rib_ts.csv'))


markers = markers[markers$p_val_adj < 0.05, ]

wt_markers = markers[markers$avg_log2FC > 0.5, ]
rib_markers = markers[markers$avg_log2FC < -0.5, ]

write.csv(wt_markers, file = file.path(TARGET_dir, 'markers_wt_sg.csv'))
write.csv(rib_markers, file = file.path(TARGET_dir, 'markers_rib_sg.csv'))

my_sample_col = data.frame("sample" = object@meta.data$experimental_condition, 
                           row.names = row.names(object@meta.data))

markers = markers[order(markers$avg_log2FC), ]
exp_dat = object@assays$RNA@data[c(rownames(wt_markers), rownames(rib_markers)), ]
exp_dat = as.matrix(exp_dat)
exp_dat = apply(exp_dat, 1, scale)
exp_dat = t(exp_dat)
colnames(exp_dat) = rownames(my_sample_col)



p = pheatmap(exp_dat, 
             annotation_col = my_sample_col, 
             show_colnames = FALSE)

avg_exp_df = matrix(nrow = nrow(exp_dat), ncol = 2)
colnames(avg_exp_df) = c('Wt', 'Rib')
rownames(avg_exp_df) = rownames(exp_dat)
avg_exp_df = as.data.frame(avg_exp_df)

exp_dat = exp_dat[, rownames(my_sample_col)]
for(temp_samp in unique(my_sample_col$sample)) {
  temp_exp_df = exp_dat[, my_sample_col$sample == temp_samp]
  
  for(gene in rownames(avg_exp_df)) {
    avg_exp_df[gene, temp_samp] = mean(temp_exp_df[gene, ])
  }
}

p = pheatmap(avg_exp_df, cluster_cols = FALSE, cluster_rows = FALSE)

# save the heatmap 
enrichR::setEnrichrSite("FlyEnrichr")

enrichment_results = enrichR::enrichr(
  genes = rownames(wt_markers), 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018", 
    'KEGG_2019'
  )
)

biological_analysis = enrichment_results$GO_Biological_Process_2018
molecular_analysis = enrichment_results$GO_Molecular_Function_2018

biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]

write.csv(biological_analysis, file = file.path(TARGET_dir, 'wt_biological_GO.csv'))
write.csv(molecular_analysis, file = file.path(TARGET_dir, "wt_molecular_function_GO.csv"))

# this is to get the rib 
enrichment_results = enrichR::enrichr(
  genes = rownames(rib_markers), 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018", 
    'KEGG_2019'
  )
)

biological_analysis = enrichment_results$GO_Biological_Process_2018
molecular_analysis = enrichment_results$GO_Molecular_Function_2018

biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
biological_analysis = biological_analysis[order(biological_analysis$Overlap, decreasing = TRUE), ]
write.csv(biological_analysis, file = file.path(TARGET_dir, 'rib_biological_GO.csv'))
write.csv(molecular_analysis, file = file.path(TARGET_dir, "rib_molecular_function_GO.csv"))

