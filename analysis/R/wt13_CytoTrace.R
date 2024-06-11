library(Matrix)
library(reticulate)
library(Seurat)
library(harmony)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "wt13_cytoTrace")
dir.create(TARGET_dir)

wt_object = readRDS(file.path('results', ANALYSIS_VERSION, 'manual_annotation_wt13/manual_celltype_object2.rds'))

withr::with_dir(TARGET_dir, {
  for(cell_type in unique(wt_object$manual_celltypes)) {
    cell_type_folder = stringr::str_replace_all(cell_type, pattern = "/", "_")
    print(cell_type_folder)
    dir.create(cell_type_folder)
    
    sub_object = subset(wt_object, subset = manual_celltypes == cell_type)
    
    query_exp = sub_object@assays$RNA@counts
    write(colnames(query_exp), file = file.path(cell_type_folder, "raw_query_colnames.txt"))
    write(rownames(query_exp), file = file.path(cell_type_folder, "raw_query_rownames.txt"))
    Matrix::writeMM(query_exp, file.path(cell_type_folder, "raw_query_exp.txt"))
    
    query_meta_tab = sub_object@meta.data
    write.csv(query_meta_tab, file = file.path(cell_type_folder, 'query_sample_tab.csv'))
  }
})

reticulate::use_condaenv("OneCC_dev")
py_run_file("Python/run_cytoTrace.py")

withr::with_dir(TARGET_dir, {

  for(cell_type in unique(wt_object$manual_celltypes)) {
    cell_type_folder = stringr::str_replace_all(cell_type, pattern = "/", "_")
    sub_object = subset(wt_object, subset = manual_celltypes == cell_type)
    
    if(nrow(sub_object@meta.data) < 100) {
      next()
    }
    print(cell_type_folder)
    
    sub_object = NormalizeData(sub_object)
    sub_object <- FindVariableFeatures(sub_object, selection.method = "vst", nfeatures = 2000)
    sub_object <- ScaleData(sub_object)
    sub_object <- RunPCA(sub_object, features = VariableFeatures(object = sub_object), npc = 20)
    
    meta_tab = read.csv(file.path(cell_type_folder, "cytotraced_sample_tab.csv"), row.names = 1)
    sub_object@meta.data$ct_pt = meta_tab$ct_pseudotime
    sub_object = RunUMAP(sub_object, dims = 1:20)
    
    sub_object %<>% harmony::RunHarmony("batch")
    sub_object = RunUMAP(sub_object, dims = 1:20, reduction = "harmony")
    
    FeaturePlot(sub_object, features = 'ct_pt')
    DimPlot(sub_object, group.by = 'batch')
    saveRDS(sub_object, file = file.path(cell_type_folder, 'sub_object.rds'))
    
    p = FeaturePlot(sub_object, features = 'ct_pt')
    ggsave(file.path(cell_type_folder, 'sub_object_ct.png'))
  }
  
})


