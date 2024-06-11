# Set up the data
wt_rep1  = readRDS(file.path("results", ANALYSIS_VERSION, "wt_rep1/object.Rdata"))
rib_rep1 = readRDS(file.path("results", ANALYSIS_VERSION, "rib_rep1/object.Rdata"))
# rib_rep2 = readRDS(file.path("results", ANALYSIS_VERSION, "rib_rep2/object.Rdata"))
wt_rep1$genotype  = "WT"
rib_rep1$genotype = "rib KO"
# rib_rep2$genotype = "rib KO"
wt_rep1$sample    = "WT rep1"
rib_rep1$sample   = "rib KO rep1"
# rib_rep2$sample   = "rib KO rep2"
object = merge(wt_rep1, rib_rep1)
# object = merge(object, rib_rep2)

#' Save joint clustering and embedding results to the current working directory
#'
saveJointEmbeddings = function(object){
  write.csv( cbind( object[["joint_cluster"]], object[["umap"]]@cell.embeddings ), 
             "umap_and_clustering.csv")
}

#' See if joint clustering and embedding results were saved already
#'
isPresentJointEmbeddings = function(){
  file.exists("umap_and_clustering.csv")
}

#' Load joint clustering and embedding results from the current working directory
#'
loadJointEmbeddings = function(object){
  object
  X = read.csv( "umap_and_clustering.csv", row.names = 1 )
  object[["umap"]] = Seurat::CreateDimReducObject(embeddings = as.matrix(X[2:3]), key = "UMAP_")
  object[["joint_cluster"]] = X[[1]]
  return(object)
}

#' Check on how the samples mix. 
#'
compareJointEmbeddings = function(object){
  object[["shortAnnotation"]][[1]] %<>% toupper
  object[["sample_and_celltype"]] = paste0( object[["sample"]][[1]], " ",  object[["shortAnnotation"]][[1]] )
  
  # Check overlap in old and new clusters
  table(new_cluster = object$joint_cluster, annotation = object$annotation, sample = object$sample) %>%
    data.frame() %>%
    dplyr::arrange(-Freq) %>%
    dplyr::arrange(new_cluster) %>%
    subset(Freq >= 10) %>%
    write.csv("New and old cluster overlap.csv")
  
  # Display old clusters on umap
  {
    # Horrible cluttered plot with sample X celltype on joint embedding
    # sparser / more readable plots follow below. 
    pdf(width = 16, height = 7, "dense_comparison.pdf")
    p = DimPlot(object, group.by = "sample_and_celltype", label = T, repel = T, reduction = "umap") + coord_fixed()
    print(p)
    dev.off()
    
    # Look for false negatives -- should be together and aren't
    pdf(width = 16, height = 7, "comparison_grouped_by_prior_labels.pdf")
    for( sa in unique(object[["shortAnnotation"]][[1]])){
      title = make.names(paste0("Prior_annotation_is_", sa))
      object[[title]] = "other cell types"
      object[[title]][object$shortAnnotation == sa,] = object[["sample"]][object$shortAnnotation == sa,]
      
      p = DimPlot(object, 
                  cols = c("other cell types" = "gray", "rib KO rep1" = "red",  "WT rep1" = "blue"),
                  group.by = title, 
                  label = T, 
                  reduction = "umap", 
                  dims = 1:2,
                  repel = T) + coord_fixed()
      print(p)
    }
    dev.off()
    # Look for false positives -- together and shouldn't be
    pdf(width = 16, height = 7, "comparison_grouped_by_new_labels.pdf")
    for( jc in unique(object[["joint_cluster"]][[1]])){
      plot_title = paste0("joint_cluster_is_", jc) 
      object[[plot_title]] = "other clusters from joint clustering"
      # The plot is cluttered and unreadable unless you remove less abundant cell types. 
      well_represented = 
        object[["sample_and_celltype"]][[1]][object$joint_cluster == jc] %>%
        table %>% 
        as.data.frame() %>% 
        subset(Freq > 50, select = ".", drop = T) 
      if(length(well_represented)<1){
        next
      }
      gets_colored = 
        object$sample_and_celltype %in% well_represented &
        object$joint_cluster == jc
      object[[plot_title]][gets_colored,] = object[["sample_and_celltype"]][gets_colored,]
      p = DimPlot(object, 
                  group.by = plot_title, 
                  cols = setNames(
                    scales::hue_pal()(length(well_represented)) %>% c("gray"),
                    well_represented %>% as.character %>% c("other clusters from joint clustering")
                  ),
                  label = T, reduction = "umap", 
                  repel = T) + coord_fixed()
      print(p)
    }
    dev.off()
  }
  # display basics on new UMAP
  MakeBasicPlots(object, plotname = "joint_embedding.pdf")
  p = DimPlot(object, group.by = "joint_cluster", reduction = "umap", label = T, repel = T) + coord_fixed()
  ggsave(width = 16, height = 7, filename = "new_clusters.pdf", plot = p)
}

#' Run a data integration or batch correction algorithm.
#'
makeJointEmbeddings = function(object, method){
  dir.create(file.path(RESULTS, method), recursive = T, showWarnings = F)
  withr::with_dir(file.path(RESULTS, method), {  
    if(!isPresentJointEmbeddings()){
      if(method=="no_batch_correction"){
        object %<>% Seurat::NormalizeData()
        object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
        object %<>% ScaleData(features = VariableFeatures(object = object), vars.to.regress = "sample")
        object %<>% RunPCA(features = VariableFeatures(object = object), npcs = 100)
        object %<>% RunUMAP(features = VariableFeatures(object = object), npcs = 100)
        object %<>% FindNeighbors()
        object %<>% FindClusters()
        object[["joint_cluster"]] = object[["seurat_clusters"]]
      } else if(method=="seurat"){
        anchors = Seurat::FindIntegrationAnchors(
          object.list = list(
            subset(object, genotype == "WT"),
            subset(object, genotype == "rib KO")
          ),
          dims = 1:100,
          anchor.features = VariableFeatures(object)
        )
        object = IntegrateData(anchorset = anchors, dims = 1:100) 
        DefaultAssay(object) <- "integrated"
        object %<>% ScaleData()
        object %<>% RunPCA(n_pcs = 100)
        object %<>% FindNeighbors()
        object %<>% FindClusters(resolution = 1)
        object[["joint_cluster"]] = object[["seurat_clusters"]]
        object %<>% RunUMAP(dims = 1:50)
      } else if(method=="harmony"){
        object %<>% Seurat::FindVariableFeatures()
        object %<>% ScaleData()
        object %<>% RunPCA(npcs = 100)
        # Set project.dim = F due to #86
        # https://github.com/immunogenomics/harmony/issues/86
        object %<>% harmony::RunHarmony(group.by.vars = c("sample"), dims.use = 1:100, assay.use = "RNA", project.dim = F)
        object %<>% FindNeighbors(reduction = "harmony")
        object %<>% FindClusters(resolution = 2, dims = 1:100, reduction = "harmony")
        object[["joint_cluster"]] = object[["seurat_clusters"]]
        object %<>% RunUMAP(dims = 1:100, reduction = "harmony")
      } else if(method=="harmony2"){
        object %<>% Seurat::FindVariableFeatures()
        object %<>% ScaleData()
        object %<>% RunPCA(npcs = 100)
        # Set project.dim = F due to #86
        # https://github.com/immunogenomics/harmony/issues/86
        object %<>% harmony::RunHarmony(group.by.vars = c("sample", "barcode_contains"),
                                        dims.use = 1:100, assay.use = "RNA", project.dim = F)
        object %<>% FindNeighbors(reduction = "harmony")
        object %<>% FindClusters(resolution = 2, dims = 1:100)
        object %<>% RunUMAP(dims = 1:100, reduction = "harmony")
        object[["joint_cluster"]] = object[["seurat_clusters"]]
      } else if(method=="harmony3"){
        object %<>% Seurat::FindVariableFeatures()
        object %<>% ScaleData()
        object %<>% RunPCA(npcs = 100)
        # Set project.dim = F due to #86
        # https://github.com/immunogenomics/harmony/issues/86
        object %<>% harmony::RunHarmony(group.by.vars = c("sample", "Phase"),
                                        dims.use = 1:100, assay.use = "RNA", project.dim = F)
        object %<>% FindNeighbors(reduction = "harmony")
        object %<>% FindClusters(resolution = 2, dims = 1:100)
        object %<>% RunUMAP(dims = 1:100, reduction = "harmony")
        object[["joint_cluster"]] = object[["seurat_clusters"]]
      } else if(method=="liger"){
        liger <- rliger::createLiger(raw.data = Seurat::SplitObject(object, split.by = "sample") %>% lapply(GetAssayData))
        liger %<>% rliger::normalize()
        liger %<>% rliger::selectGenes(var.thresh = 0.1)
        liger %<>% rliger::scaleNotCenter()
        liger %<>% rliger::optimizeALS(k=30, use.unshared = TRUE, max_iters =30,thresh=1e-10)
        liger %<>% rliger::quantile_norm()
        liger %<>% rliger::louvainCluster()
        liger %<>% rliger::runUMAP()
        object[["joint_cluster"]] = liger@clusters[Cells(object)] 
        object[["umap"]] = SeuratObject::CreateDimReducObject(embeddings = liger@tsne.coords[Cells(object), ], key = "umap") 
      } else {
        stop("No method by that name.\n")
      }
      saveJointEmbeddings(object)
    } else {
      object = loadJointEmbeddings(object)
    }
    compareJointEmbeddings(object)
  })
  object
}


object %<>% makeJointEmbeddings("no_batch_correction")
object %<>% makeJointEmbeddings("seurat")
object %<>% makeJointEmbeddings("harmony")
object %<>% makeJointEmbeddings("harmony2")
object %<>% makeJointEmbeddings("harmony3")
object %<>% makeJointEmbeddings("liger")

# How consistent is the difference between clusters?
# 
{
  X = list()
  X[["all"]] = FindMarkers(object, 
                           ident.1 = "WT rep1",
                           ident.2 = "rib KO rep1",  
                           group.by = "sample",
                           min.pct = 0, 
                           min.cells.group = 0, 
                           logfc.threshold = 0, 
                           test.use = "bimod")
  for( sa in unique(object[["shortAnnotation"]][[1]])){
    try({
      X[[sa]] = FindMarkers(object[, object$shortAnnotation == sa], 
                          ident.1 = "WT rep1",
                          ident.2 = "rib KO rep1", 
                          group.by = "sample" ,
                          min.pct = 0, 
                          min.cells.group = 0, 
                          logfc.threshold = 0, 
                          test.use = "bimod")
    })
  }
  dir.create(file.path(RESULTS, "differential_expression"), recursive = F, showWarnings = T)
  mapply(
    function(X, name) { write.csv(X, file.path(RESULTS, "differential_expression", paste0(name, ".csv"))) },
    X, 
    names(X)
  )
  genes = rownames(X[["all"]])
  {
    pdf(file.path(RESULTS, "differential_expression", "fold_change_per_cluster.pdf"))
    for( sa in names(X)){
      plot(X[["all"]][genes, "avg_log2FC"],
           X[[sa]][genes, "avg_log2FC"], 
           pch = ".",
           xlab = paste0("Fold change cluster ", sa),
           ylab = "Fold change over all clusters",
           main = sa) 
      abline(a=0, b=1, col = "red")
    }
    dev.off()
  }
}
differential_abundance = table(object$shortAnnotation, object$sample)
differential_abundance %>%
  as.matrix.data.frame() %>% 
  set_rownames(rownames(differential_abundance)) %>%
  set_colnames(colnames(differential_abundance)) %>%
  write.csv(file.path(RESULTS, "differential_abundance.csv"))
