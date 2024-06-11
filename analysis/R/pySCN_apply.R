library(singleCellNet)
library(reticulate)
library(Matrix)
library(Seurat)

RESULTS = file.path("results", ANALYSIS_VERSION, 'SCN_application')
dir.create(RESULTS)

withr::with_dir(
  file.path(RESULTS), 
  {
    for(i in rownames(metadata)){ 
      sample_id = metadata[i, 'sample']
      print(sample_id)
      seurat_obj = readRDS(file.path("..", sample_id, 'object.Rdata'))
      dir.create(sample_id)
      
      query_exp = seurat_obj@assays$RNA@counts
      write(colnames(query_exp), file = file.path(sample_id, "raw_query_colnames.txt"))
      write(rownames(query_exp), file = file.path(sample_id, "raw_query_rownames.txt"))
      Matrix::writeMM(query_exp, file.path(sample_id, "raw_query_exp.txt"))
    }
  }
)

# run the python script 
# python run_pySCN.py V18 
withr::with_dir(
  file.path(RESULTS), 
  {
    for(i in rownames(metadata)){ 
      sample_id = metadata[i, 'sample']
      print(sample_id)
      seurat_obj = readRDS(file.path("..", sample_id, 'object.Rdata'))
      SCN_class_tab = read.csv(paste0(sample_id, "/SCN_classification.csv"), row.names = 1)
      seurat_obj@meta.data = cbind(seurat_obj@meta.data, SCN_class_tab)
      
      for(temp_col in colnames(SCN_class_tab)) { 
        if(temp_col == 'SCN_class') { 
          p = DimPlot(seurat_obj, group.by = temp_col)
        } else { 
          p = FeaturePlot(seurat_obj, features = temp_col)
        }
        ggsave(file.path(sample_id, paste0(temp_col,'_SCN_results.png')), plot = p, height = 6, width = 8)
      }
    }
  }
)
