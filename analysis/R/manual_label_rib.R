library(Seurat)
library(harmony)
library(magrittr)
library(ggplot2)
library(stringr)
library(reticulate)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "manual_annotation_rib")
dir.create(TARGET_dir)

train_obj = readRDS(file.path("results", ANALYSIS_VERSION, 'manual_annotation_wt13/manual_celltype_object2.rds'))
withr::with_dir(TARGET_dir, {
  dir.create('training_data')
  raw_expTrain = train_obj@assays$RNA@counts
  raw_stTrain = train_obj@meta.data
  write(colnames(raw_expTrain), file = "training_data/raw_train_colnames.txt")
  write(rownames(raw_expTrain), file = "training_data/raw_train_rownames.txt")
  Matrix::writeMM(raw_expTrain, "training_data/raw_train_exp.txt")
  write.csv(raw_stTrain, file = 'training_data/raw_meta_tab.csv')
})

query_obj = readRDS(file.path("results", ANALYSIS_VERSION, "rib_integrated/BDGP_automated_annotation_object.rds"))
withr::with_dir(TARGET_dir, {
  dir.create('query_data')
  raw_expQuery = query_obj@assays$RNA@counts
  raw_stQuery = query_obj@meta.data
  write(colnames(raw_expQuery), file = "query_data/raw_query_colnames.txt")
  write(rownames(raw_expQuery), file = "query_data/raw_query_rownames.txt")
  Matrix::writeMM(raw_expQuery, "query_data/raw_query_exp.txt")
  write.csv(raw_stQuery, file = 'query_data/raw_meta_tab.csv')
})

# recluster the results 
query_obj %<>% FindClusters(resolution = 1)
query_obj %<>% FindClusters(resolution = 2)

DimPlot(query_obj, label = TRUE)

withr::with_dir(TARGET_dir, {
  saveRDS(query_obj, file = 'new_clustering_obj.rds')
  p = DimPlot(query_obj, label = TRUE, group.by = "RNA_snn_res.1")
  ggsave(filename = file.path('cluster_UMAP.png'), plot = p, width = 12, height = 9)
})

withr::with_dir(
  file.path(TARGET_dir), 
  {
    if('top_marker_genes_1-20.csv' %in% list.files() == FALSE) { 
      differentially_expressed = FindAllMarkers(query_obj, test.use = "bimod")
      write.csv(differentially_expressed, "marker_genes.csv")
      differentially_expressed %>%
        dplyr::group_by(cluster) %>%
        dplyr::top_n(n=20, wt = avg_log2FC) %>%
        dplyr::arrange(cluster) %>%
        write.csv("top_marker_genes_1-20.csv")      
    } else { 
      differentially_expressed = read.csv("top_marker_genes_1-20.csv", row.names = 1)
    }
  }
)

################################
# this is to work on the marker genes previously defined 
marker_genes = read.csv(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13", 'manualMarkerGenes.csv'))
marker_genes_list = list()

for(unique_ct in unique(marker_genes$annotation)) {
  sub_marker_genes = marker_genes[marker_genes$annotation == unique_ct, ]
  marker_genes_list[[trimws(unique_ct)]] = intersect(as.vector(sub_marker_genes$marker_genes), rownames(query_obj@assays$RNA@data))
  query_obj = AddModuleScore(
    object = query_obj,
    features = list(marker_genes_list[[trimws(unique_ct)]]),
    name = paste0(stringr::str_replace_all(trimws(unique_ct), pattern = " ", replacement = "_"), '_ModuleScore')
  )
}

module_scores_names = colnames(query_obj@meta.data)[grep("ModuleScore", colnames(query_obj@meta.data))]

withr::with_dir(TARGET_dir, {
  dir.create('moduleScores_UMAP')
  dir.create('moduleScores_Vln')
  
  for(module_name in module_scores_names) {
    p = FeaturePlot(query_obj, features = module_name)
    ggsave(paste0('moduleScores_UMAP/', module_name, ".png"), plot = p, width = 8, height = 6)
    
    p = VlnPlot(query_obj, features = module_name, pt.size = 0, group.by = 'RNA_snn_res.1')
    ggsave(paste0('moduleScores_Vln/', module_name, ".png"), plot = p, width = 8, height = 6)
  }
})

####################
names(table(raw_stTrain$manual_celltypes))
query_obj@meta.data$manual_celltypes = NA
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(0), 'manual_celltypes'] = "something"
DimPlot(query_obj, label = TRUE, group.by = 'manual_celltypes')

names(table(raw_stTrain$manual_celltypes))
query_obj@meta.data$manual_celltypes = NA
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(19, 21), 'manual_celltypes'] = "Fat Body"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(4, 5, 10), 'manual_celltypes'] = "Somatic Muscle"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(18), 'manual_celltypes'] = "Dorsal Vessel (heart)"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(14), 'manual_celltypes'] = "Visceral Mesoderm"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(8), 'manual_celltypes'] = "Muscle System"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(31), 'manual_celltypes'] = "Longitudinal/Caudal Visceral Mesoderm"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(24), 'manual_celltypes'] = "Plasmatocytes"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(9, 23), 'manual_celltypes'] = "Glia"

# these are crystal cell genes 
gene = 'PPO2'
p = FeaturePlot(query_obj, features = gene)
ggsave(file.path(TARGET_dir, paste0(gene, '_exp.png')), plot = p, width = 5, height = 4)

gene = 'PPO1'
p = FeaturePlot(query_obj, features = gene)
ggsave(file.path(TARGET_dir, paste0(gene, '_exp.png')), plot = p, width = 5, height = 4)
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(26), 'manual_celltypes'] = "Crystal Cells"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(28), 'manual_celltypes'] = "Malpighian Tubules"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(20, 7), 'manual_celltypes'] = "Gut Endoderm"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(29), 'manual_celltypes'] = "Germ Cell"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(16, 3), 'manual_celltypes'] = "Early Cells"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(22), 'manual_celltypes'] = "Yolk"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(2, 13, 1, 11), 'manual_celltypes'] = "CNS"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(15), 'manual_celltypes'] = "Head Sensory System"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(27), 'manual_celltypes'] = "Amnioserosa"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(30), 'manual_celltypes'] = "Glia Midline" 
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(6), 'manual_celltypes'] = "Optic Lobe"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(17), 'manual_celltypes'] = "Dorsal Epidermis"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(0), 'manual_celltypes'] = "Ventral Epidermis"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(12), 'manual_celltypes'] = "Trachea"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(32), 'manual_celltypes'] = "Garland Cells"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.1 %in% c(25), 'manual_celltypes'] = "Unknown"

withr::with_dir(TARGET_dir, {
  saveRDS(query_obj, file = 'manual_celltyping_1_obj.rds')
  p = DimPlot(query_obj, label = TRUE, group.by = "manual_celltypes")
  ggsave(filename = file.path('manual_celltypes1.png'), plot = p, width = 12, height = 9)
})

# then we proceed with second round of cell typing 
names(table(raw_stTrain$manual_celltypes))
DimPlot(query_obj, group.by = 'RNA_snn_res.2', label = TRUE)
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.2 %in% c(28), 'manual_celltypes'] = "Lateral Sensory Neurons"
query_obj@meta.data[query_obj@meta.data$RNA_snn_res.2 %in% c(35), 'manual_celltypes'] = "Hindgut"

withr::with_dir(TARGET_dir, {
  saveRDS(query_obj, file = 'manual_celltyping_2_obj.rds')
  p = DimPlot(query_obj, label = TRUE, group.by = "manual_celltypes")
  ggsave(filename = file.path('manual_celltypes2.png'), plot = p, width = 15, height = 9)
})

# find the rare SG population by subclustering
source('R/manual_find_rib_SG.R')

find_rib_obj = readRDS(file.path("results", ANALYSIS_VERSION, "manual_find_rib_SG", 'rib_find_SG_object.rds'))
SG_cells = rownames(find_rib_obj@meta.data[find_rib_obj@meta.data$finer_annotation == 'salivary gland', ])
query_obj@meta.data[SG_cells, 'manual_celltypes'] = 'Salivary Gland'

withr::with_dir(TARGET_dir, {
  saveRDS(query_obj, file = 'manual_celltyping_3_obj.rds')
  p = DimPlot(query_obj, label = TRUE, group.by = "manual_celltypes")
  ggsave(filename = file.path('manual_celltypes3.png'), plot = p, width = 15, height = 9)
})

query_obj = readRDS(file.path(TARGET_dir, 'manual_celltyping_3_obj.rds'))
query_obj@meta.data[query_obj@meta.data$manual_celltypes == 'Esophagus', 'manual_celltypes'] = 'Epidermis'
query_obj@meta.data[query_obj@meta.data$manual_celltypes == 'Ventral Epidermis', 'manual_celltypes'] = 'Epidermis'
query_obj@meta.data[query_obj@meta.data$manual_celltypes == 'Hypopharynx', 'manual_celltypes'] = 'Epidermis'
query_obj@meta.data[query_obj@meta.data$manual_celltypes == 'Dorsal Epidermis', 'manual_celltypes'] = 'Epidermis'
query_obj@meta.data[query_obj@meta.data$manual_celltypes == 'Epipharynx', 'manual_celltypes'] = 'Epidermis'

p = DimPlot(query_obj, group.by = 'manual_celltypes', label = TRUE)
ggsave(filename = file.path(TARGET_dir, 'manual_celltypes4.png'), plot = p, width = 12, height = 8)
saveRDS(query_obj, file.path(TARGET_dir, 'manual_celltyping_4_obj.rds'))

####################

query_obj = readRDS(file.path(TARGET_dir, 'manual_celltyping_4_obj.rds'))
# this is to run the pySCN 
py_run_file(file.path(TARGET_dir, 'run_scn_python.py'))

rib_class_results = read.csv(file.path(TARGET_dir, 'pySCN_results/rib_SCN_classification.csv'), row.names = 1)
rib_class_results = rib_class_results[rownames(query_obj@meta.data), ]
query_obj@meta.data = cbind(query_obj@meta.data, rib_class_results)

p = DimPlot(query_obj, group.by = 'SCN_class', label = TRUE)
ggsave(filename = file.path(TARGET_dir, 'SCN_class.png'), plot = p, width = 18, height = 9)

p = FeaturePlot(query_obj, features = 'Salivary.Gland')
ggsave(filename = file.path(TARGET_dir, 'SG_class.png'), plot = p, width = 10, height = 9)

p = FeaturePlot(query_obj, features = 'Salivary_Gland_ModuleScore1')
ggsave(filename = file.path(TARGET_dir, 'SG_gene_module.png'), plot = p, width = 10, height = 9)


meta_tab = query_obj@meta.data

plot_df = meta_tab %>%
  group_by(manual_celltypes, SCN_class) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

ggplot(plot_df, aes(x = manual_celltypes, y = freq, fill = SCN_class)) + 
  geom_bar(stat = "identity")

withr::with_dir(TARGET_dir, {
  dir.create('SCN_proportion')
  for(cluster in unique(plot_df$manual_celltypes)) {
    temp_df = plot_df[plot_df$manual_celltypes == cluster, ]
    
    p = ggplot(temp_df, aes(x = reorder(SCN_class, -freq), y = freq, fill = SCN_class)) + 
      geom_bar(stat = "identity") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1), legend.position="none") + 
      xlab('SCN_class') + 
      ylim(c(0, 1)) + 
      ggtitle(paste0('cluster ', cluster))
    
    ggsave(filename = paste0('SCN_proportion/', cluster, '.png'), plot = p, width = 8, height = 6)
  }
})












