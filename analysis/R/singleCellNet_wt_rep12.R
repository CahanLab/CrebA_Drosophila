library(Seurat)
library(stringr)
library(singleCellNet)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "wt_trained_SCN_rep12")
dir.create(TARGET_dir)

# this is specifically running for V18 
wt_1_seurat = readRDS("results/v18/wt_rep1/object.Rdata")
wt_2_seurat = readRDS("results/v18/wt_rep2/object.Rdata")

wt_1_labels = read.csv("results/v18/Bianca_wt_cluster/wt1_clusters.csv")
colnames(wt_1_labels) = c('cluster_number', 'cluster_identification')

wt_2_labels = read.csv("results/v18/Bianca_wt_cluster/wt2_clusters.csv")
colnames(wt_2_labels) = c('cluster_number', 'cluster_identification')

# assigning the cell-types 
assign_ct <- function(seurat_obj, annotate_st) { 
  seurat_obj@meta.data$Bianca_CT = NA
  for(cluster in unique(annotate_st$cluster_number)) { 
    seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters == str_trim(cluster), 'Bianca_CT'] = annotate_st[annotate_st$cluster_number == cluster, 'cluster_identification']
  }
  return(seurat_obj)
}

wt_1_seurat = assign_ct(wt_1_seurat, wt_1_labels)
wt_1_seurat$batch = 'wt_1'

wt_2_seurat = assign_ct(wt_2_seurat, wt_2_labels)
wt_2_seurat$batch = 'wt_2'

combined_seurat = merge(x = wt_1_seurat, y = wt_2_seurat, project = 'batch')
common_ct_list = intersect(wt_1_seurat$Bianca_CT, wt_2_seurat$Bianca_CT)

subset_seurat <- subset(x = combined_seurat, subset = Bianca_CT %in% c(common_ct_list, 'Epidermis'))

raw_expTrain = subset_seurat@assays$RNA@counts
raw_stTrain = subset_seurat@meta.data

library(Matrix)
withr::with_dir(
  file.path(TARGET_dir), 
  { 
    
    write(colnames(raw_expTrain), file = "raw_train_colnames.txt")
    write(rownames(raw_expTrain), file = "raw_train_rownames.txt")
    Matrix::writeMM(raw_expTrain, "raw_train_exp.txt")
    
    write.csv(raw_stTrain, file = 'raw_meta_tab.csv')
    
  }
)

# run the python script to train the classifier 
query_obj = readRDS(file.path('results', ANALYSIS_VERSION, "rib_integrated/object.rds"))
query_exp = query_obj@assays$RNA@counts

withr::with_dir(
  file.path(TARGET_dir), 
  {
    write(colnames(query_exp), file = file.path("raw_query_colnames.txt"))
    write(rownames(query_exp), file = file.path("raw_query_rownames.txt"))
    Matrix::writeMM(query_exp, file.path("raw_query_exp.txt"))
  }
)

# run the python script to apply the classifier 
withr::with_dir(
  file.path(TARGET_dir), 
  {
    SCN_class_tab = read.csv(paste0("rib_SCN_classification.csv"), row.names = 1)
    query_obj@meta.data = cbind(query_obj@meta.data, SCN_class_tab)
    dir.create('rib_integrated_results')
    
    for(temp_col in colnames(SCN_class_tab)) { 
      if(temp_col == 'SCN_class') { 
        p = DimPlot(query_obj, group.by = temp_col)
      } else { 
        p = FeaturePlot(query_obj, features = temp_col)
      }
      ggsave(file.path('rib_integrated_results', paste0(temp_col,'_SCN_results.png')), plot = p, height = 6, width = 8)
    }
    saveRDS(query_obj, file = 'classified_query_obj.rds')
  }
)

# get the sg classified cells from query 
sg_query_obj = subset(query_obj, subset = salivary.gland > 0.4)
sg_query_obj = subset(sg_query_obj, subset = SCN_class == 'salivary gland')

sg_query_obj@meta.data$experimental_condition = 'Rib'

# wt obj 
wt_obj = readRDS(file.path("results", ANALYSIS_VERSION, 'wt_integrated/object.rds'))
sg_wt_obj = subset(wt_obj, subset = Bianca_CT == unique(wt_obj@meta.data$Bianca_CT)[25]) # subset out the salivary gland
sg_wt_obj@meta.data$experimental_condition = 'Wt'

combined_obj = merge(sg_query_obj, sg_wt_obj)
sg_raw_object = Seurat::CreateSeuratObject(combined_obj@assays$RNA@counts, project = 'sg_cells')
sg_raw_object@meta.data$experimental_condition = combined_obj@meta.data$experimental_condition

sg_norm_object = Seurat::NormalizeData(sg_raw_object)

markers <- FindMarkers(object = sg_norm_object, ident.1 = "Wt", group.by = 'experimental_condition', test.use = 'wilcox') 
#markers = markers[markers$avg_log2FC < 0, ]

write.csv(markers, file = file.path(TARGET_dir, 'DE_genes_wt_rib_sg.csv'))
markers = markers[markers$p_val_adj < 0.05, ]

library(enrichR)
enrichR::setEnrichrSite("FlyEnrichr")

enrichment_results = enrichR::enrichr(
  genes = rownames(markers[markers$avg_log2FC < 0, ]), 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018"
  )
)

biological_analysis = enrichment_results$GO_Biological_Process_2018
molecular_analysis = enrichment_results$GO_Molecular_Function_2018

write.csv(biological_analysis, file = file.path(TARGET_dir, 'biological_GO.csv'))
write.csv(molecular_analysis, file = file.path(TARGET_dir, "molecular_function_GO.csv"))
