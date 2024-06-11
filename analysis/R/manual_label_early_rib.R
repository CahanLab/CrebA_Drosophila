TARGET_dir = file.path("results", ANALYSIS_VERSION, "manual_annotation_rib_early")
dir.create(TARGET_dir, recursive = TRUE)
object = readRDS(file.path("results", ANALYSIS_VERSION, "rib_early_integrated/BDGP_automated_annotation_object.rds"))
manual_annotation = read.csv(file.path("accessory_data/manual_annotation", ANALYSIS_VERSION, "rib_early/manual_annotations.csv"))

object@meta.data$manual_celltypes = NA
for(cluster in unique(object@meta.data$seurat_clusters)) {
  object@meta.data[object@meta.data$seurat_clusters == cluster, 'manual_celltypes'] = manual_annotation[manual_annotation$Cluster_id == cluster, 'Manual_Annotation']
}

DimPlot(object, group.by = 'manual_celltypes', label = TRUE)
ggsave(file.path(TARGET_dir, "manual_label.png"), width = 12, height = 8)

DimPlot(object, group.by = 'seurat_clusters', label = TRUE)

saveRDS(object, file = file.path(TARGET_dir, 'manual_celltype_object.rds'))
