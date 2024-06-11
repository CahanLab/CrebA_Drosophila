#' Re-analyze a subset of a Seurat object.
#'
runSubClustering = function(object, cluster){
  # To record results
  cluster_string = paste0(cluster, collapse = "+")
  dir.create(file.path(RESULTS, "subcluster", cluster_string), recursive = T, showWarnings = F)
  
  # Get subset
  cluster0 = subset(object, seurat_clusters %in% cluster)
  
  # Basic analysis, determining number of PC's with M-P upper bound
  cluster0 %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  cluster0 %<>% ScaleData(features = VariableFeatures(object = cluster0))
  cluster0 %<>% RunPCA(features = VariableFeatures(object = cluster0), npcs = 100)
  pc_cutoff = RMTstat::qmp( 1, ndf=length(Cells(cluster0)), pdim=length(VariableFeatures(cluster0)), var=1)
  singular_values = slot(Reductions(cluster0, slot = "pca"), "stdev")
  is_significant = singular_values^2 > pc_cutoff
  num_pc = sum(is_significant)
  cat( "Num significant PC's: ", num_pc, file = file.path(RESULTS, "subcluster", cluster_string, "pca_test.txt" ) )
  cluster0 %<>% FindNeighbors()
  cluster0 %<>% FindClusters()
  
  # plot with previous umap
  withr::with_dir(file.path(RESULTS, "subcluster", cluster_string), 
                  MakeBasicPlots(cluster0,
                                 reduction = "umap", 
                                 plotname = "subset.pdf"))
  
  # plot with new umap
  cluster0 %<>% RunUMAP(dims = 1:20)
  withr::with_dir(file.path(RESULTS, "subcluster", cluster_string), 
                  MakeBasicPlots(cluster0,
                                 reduction = "umap", 
                                 plotname = "reanalysis.pdf"))
  
  write.csv(FindAllMarkers(cluster0), file.path(RESULTS, "subcluster", cluster_string, "marker_genes.csv"))
}

runSubClustering(object, object[["seurat_clusters"]][[1]] %>% unique %>% setdiff("0") %>% sort)
runSubClustering(object, 2)
runSubClustering(object, 0)
runSubClustering(object, 4)
runSubClustering(object, 5)
runSubClustering(object, 6)
runSubClustering(object, 16)
