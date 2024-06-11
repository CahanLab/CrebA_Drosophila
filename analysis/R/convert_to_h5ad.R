library(reticulate)
library(sceasy)

use_condaenv('OneCC_dev')
loompy <- reticulate::import('loompy')

TARGET_dirs = c(file.path("results", ANALYSIS_VERSION, "parse_wt_epidermis"), 
                file.path("results", ANALYSIS_VERSION, "parse_wt_CNS"), 
                file.path("results", ANALYSIS_VERSION, "parse_wt_muscle"))

for(TARGET_dir in TARGET_dirs) { 
  if(grepl('parse', TARGET_dir)) {
    object = readRDS(file.path(TARGET_dir, 'finer_annotation_object.rds'))
    sceasy::convertFormat(object, from="seurat", to="anndata",
                          outFile=file.path(TARGET_dir, 'finer_annotation_object.h5ad'))
  }
  else {
    object = readRDS(file.path(TARGET_dir, 'BDGP_automated_annotation_object.rds'))
    sceasy::convertFormat(object, from="seurat", to="anndata",
                          outFile=file.path(TARGET_dir, 'BDGP_automated_annotation_object.h5ad'))
  }
}


# the manual wt and rib 
object = readRDS("results/v18/manual_annotation_wt13/manual_celltype_object2.rds")
sceasy::convertFormat(object, from="seurat", to="anndata",
                      outFile=file.path('results/v18/manual_annotation_wt13/', 'for_cellxgenes_object.h5ad'))

object = readRDS("results/v18/manual_annotation_rib/manual_celltyping_4_obj.rds")
sceasy::convertFormat(object, from="seurat", to="anndata",
                      outFile=file.path("results/v18/manual_annotation_rib/", 'for_cellxgenes_object.h5ad'))

object = readRDS("results/v18/rib_early_integrated/BDGP_automated_annotation_object.rds")
sceasy::convertFormat(object, from="seurat", to="anndata",
                      outFile=file.path('results/v18/rib_early_integrated/', 'for_cellxgenes_object.h5ad'))

object = readRDS("results/v18/crebA_early_integrated/BDGP_automated_annotation_object.rds")
sceasy::convertFormat(object, from="seurat", to="anndata",
                      outFile=file.path('results/v18/crebA_early_integrated/', 'for_cellxgenes_object.h5ad'))
