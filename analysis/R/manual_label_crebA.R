TARGET_dir = file.path("results", ANALYSIS_VERSION, "manual_annotation_crebA")
dir.create(TARGET_dir, recursive = TRUE)
object = readRDS(file.path("results", ANALYSIS_VERSION, "crebA_rep3/BDGP_automated_annotation_object.rds"))
manual_annotation = read.csv(file.path("accessory_data/manual_annotation", ANALYSIS_VERSION, "crebA/manual_annotations.csv"))

object@meta.data$manual_celltypes = NA
for(cluster in unique(object@meta.data$seurat_clusters)) {
  object@meta.data[object@meta.data$seurat_clusters == cluster, 'manual_celltypes'] = manual_annotation[manual_annotation$Cluster_id == cluster, 'Manual_Annotation']
}

DimPlot(object, group.by = 'seurat_clusters', label = TRUE)
DimPlot(object, group.by = 'manual_celltypes', label = TRUE)
saveRDS(object, file = file.path(TARGET_dir, 'manual_celltype_object.rds'))
