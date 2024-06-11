library(Seurat)
library(stringr)
library(singleCellNet)
library(Matrix)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "wt_trained_SCN_rep13")
dir.create(TARGET_dir)

wt_seurat = readRDS(file.path("results", ANALYSIS_VERSION, "wt13_integrated", "BDGP_automated_annotation_object.rds"))
raw_expTrain = wt_seurat@assays$RNA@counts
raw_stTrain = wt_seurat@meta.data

withr::with_dir(
  file.path(TARGET_dir), 
  { 
    write(colnames(raw_expTrain), file = "raw_train_colnames.txt")
    write(rownames(raw_expTrain), file = "raw_train_rownames.txt")
    Matrix::writeMM(raw_expTrain, "raw_train_exp.txt")
    write.csv(raw_stTrain, file = 'raw_meta_tab.csv')
  }
)

# TODO run the training python script in the folder to train the classifier 
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

# TODO run the classification python script in the folder 
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
    
    # get the sg classified cells from query 
    sg_query_obj = subset(query_obj, subset = salivary.gland > 0.2)
    sg_query_obj = subset(sg_query_obj, subset = SCN_class == 'salivary gland')
    sg_query_obj@meta.data$experimental_condition = 'Rib'
    
    query_obj@meta.data$salivary_gland_call = 'other'
    query_obj@meta.data[rownames(sg_query_obj@meta.data), 'salivary_gland_call'] = 'salivary gland'
    
    saveRDS(query_obj, file = 'classified_query_obj.rds')
  }
)


# wt obj 
wt_obj = readRDS(file.path("results", ANALYSIS_VERSION, 'wt13_integrated/BDGP_automated_annotation_object.rds'))
sg_wt_obj = subset(wt_obj, subset = Integrated_tentativeCellType == 'salivary gland') # subset out the salivary gland
sg_wt_obj@meta.data$experimental_condition = 'Wt'

combined_obj = merge(sg_query_obj, sg_wt_obj)
sg_raw_object = Seurat::CreateSeuratObject(combined_obj@assays$RNA@counts, project = 'sg_cells')
sg_raw_object@meta.data$experimental_condition = combined_obj@meta.data$experimental_condition
sg_raw_object@meta.data$batch = paste0(combined_obj@meta.data$experimental_condition, "_", combined_obj@meta.data$batch)

sg_norm_object = Seurat::NormalizeData(sg_raw_object)

sg_norm_object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
sg_norm_object %<>% ScaleData(features = VariableFeatures(object = sg_norm_object))
sg_norm_object %<>% RunPCA(features = VariableFeatures(object = sg_norm_object), npcs = 100)
pc_cutoff = RMTstat::qmp( 1, ndf=length(Cells(sg_norm_object)), pdim=length(VariableFeatures(sg_norm_object)), var=1)
singular_values = slot(Reductions(sg_norm_object, slot = "pca"), "stdev")
is_significant = singular_values^2 > pc_cutoff
num_pc = sum(is_significant)

sg_norm_object %<>% RunPCA(features = VariableFeatures(object = sg_norm_object), npcs = num_pc)
sg_norm_object %<>% FindNeighbors(dims = 1:num_pc)
sg_norm_object %<>% RunUMAP(dim = 1:num_pc)
sg_norm_object %<>% FindClusters()

DimPlot(sg_norm_object)
DimPlot(sg_norm_object, group.by = 'experimental_condition')
DimPlot(sg_norm_object, group.by = 'batch')

sg_norm_object %<>% harmony::RunHarmony("batch")
sg_norm_object %<>% FindNeighbors(reduction = "harmony", dims = 1:num_pc)
sg_norm_object %<>% RunUMAP(reduction = "harmony", dim = 1:num_pc)
sg_norm_object %<>% FindClusters()

DimPlot(sg_norm_object)
DimPlot(sg_norm_object, group.by = 'experimental_condition')
DimPlot(sg_norm_object, group.by = 'batch')
FeaturePlot(sg_norm_object, features = 'nFeature_RNA')

markers <- FindMarkers(object = sg_norm_object, ident.1 = "Wt", group.by = 'experimental_condition', test.use = 'wilcox') 
#markers = markers[markers$avg_log2FC < 0, ]

write.csv(markers, file = file.path(TARGET_dir, 'DE_genes_wt_rib_sg.csv'))
markers = markers[markers$p_val_adj < 0.05, ]

library(enrichR)
enrichR::setEnrichrSite("FlyEnrichr")

enrichment_results = enrichR::enrichr(
  genes = rownames(markers[markers$avg_log2FC > 0, ]), 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018", 
    'KEGG_2019'
  )
)

biological_analysis = enrichment_results$GO_Biological_Process_2018
biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]

molecular_analysis = enrichment_results$GO_Molecular_Function_2018

KEGG_analysis = enrichment_results$KEGG_2019

write.csv(biological_analysis, file = file.path(TARGET_dir, 'biological_GO.csv'))
write.csv(molecular_analysis, file = file.path(TARGET_dir, "molecular_function_GO.csv"))





