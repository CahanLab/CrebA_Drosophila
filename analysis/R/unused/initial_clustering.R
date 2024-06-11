
withr::with_dir(file.path("results", ANALYSIS_VERSION), {
  cellranger39k = readRDS("initial_exploration.Rdata")
})

for( barcode_set in c( "both", "neither", "cells", "nuclei", "all" ) ){
  dir.create(file.path("results", ANALYSIS_VERSION, barcode_set), recursive = T)
  categories_to_take = switch(
    barcode_set, 
    all = unique(cellranger39k$barcode_contains), 
    both = c("nucleus", "cell"), 
    cells = "cell", 
    nuclei = "nucleus",
    neither = "unknown"
  )
  current_barcode_set = cellranger39k[,cellranger39k$barcode_contains %in% categories_to_take ]
  gc()
  
  # Basic analysis
  current_barcode_set %<>% Seurat::NormalizeData()
  current_barcode_set %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
  current_barcode_set %<>% ScaleData(features = VariableFeatures(object = current_barcode_set))
  current_barcode_set %<>% RunPCA(features = VariableFeatures(object = current_barcode_set), npcs = 100)
  current_barcode_set %<>% FindNeighbors()
  current_barcode_set %<>% FindClusters()
  current_barcode_set %<>% RunUMAP(dims = 1:100)
  
  # Save results
  withr::with_dir(file.path("results", ANALYSIS_VERSION, barcode_set), {
    saveRDS(current_barcode_set, "object.Rdata")
  })
  
  # Show UMAP, clusters, per-cell depth, doublet score, and various QC metrics
  withr::with_dir(file.path("results", ANALYSIS_VERSION, barcode_set), {
    DimPlot(current_barcode_set, group.by = "barcode_contains", label = T) + 
      coord_fixed() + 
      ggtitle("Droplet contents")
    ggsave("Droplet contents.pdf", width = 6, height = 6)
    
    pdf("umap_qc.pdf", width = 10, height = 8)
    DimPlot(current_barcode_set, group.by = "Phase", label = T, cols = COLORSCALES$phase) %>% print
    qc_features = c(
      "seurat_clusters",
      "DoubletFinder_call",
      "barcode_contains",
      "log10_nCount_RNA", 
      "nCount_RNA", 
      "nFeature_RNA",
      "fraction_unspliced",
      "total_transcripts_unspliced",
      "total_transcripts_spliced",
      "ribosomal_rna_norm_expression",
      "ribosomal_protein_norm_expression",
      "mitochondrial_transcript_norm_expression",
      "DoubletFinder_score", 
      "sage"
    )
    for( f in qc_features ){
      try({
        if(is.numeric(current_barcode_set[[f]][[1]])){
          FeaturePlot(current_barcode_set, features = f, coord.fixed = T) %>% print
        } else {
          DimPlot(    current_barcode_set, group.by = f, label = T) %>% print
        }
      }, silent = T)
    }
    dev.off()
  })
}  