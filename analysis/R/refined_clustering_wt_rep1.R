#' Run two analyses: one typical workflow (gene selection, PCA, clustering, UMAP)
#' and another with integration of cells, nuclei, and unknowns.  
#' 
withr::with_dir(RESULTS, {
  object = readRDS("object.Rdata")
})

old_cluster = read.csv("results/v17/wt_rep1/before_cleaning/umap_and_clustering.csv", row.names = 1)
old_cluster = old_cluster[colnames(object), ]
object@meta.data$old_cluster = old_cluster$seurat_clusters

# Redo analysis after cleaning out empties and doublets
object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
object %<>% ScaleData(features = VariableFeatures(object = object))
object %<>% RunPCA(features = VariableFeatures(object = object), npcs = 100)
pc_cutoff = RMTstat::qmp( 1, ndf=length(Cells(object)), pdim=length(VariableFeatures(object)), var=1)
singular_values = slot(Reductions(object, slot = "pca"), "stdev")
is_significant = singular_values^2 > pc_cutoff
num_pc = sum(is_significant)
object %<>% FindNeighbors(reduction = 'pca', dims = 1:num_pc) # might have to decrease the number of neighbors to find rare populations 
object %<>% FindClusters()
object %<>% RunUMAP(reduction = "pca", dims = 1:num_pc, n.components = 3L)
object[["umap3D"]] = Reductions(object, "umap")
object %<>% RunUMAP(dims = 1:num_pc)

# Show UMAP, clusters, per-cell depth, doublet score, and various QC metrics
dir.create(RESULTS)
withr::with_dir(RESULTS, {
  MakeBasicPlots(object, reduction = "umap", plotname = "not_integrated.pdf")
  save3DUMAPCategorical(object, umapkey = "umap3D", group.by = "seurat_clusters")
  saveRDS(object, "object.Rdata")  
})
