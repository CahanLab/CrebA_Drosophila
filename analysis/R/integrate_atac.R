# When we look for ribbon targets and possible ribbon cofactors, ATAC-seq could be useful.
# But it's hard to do data integration and get consistent annotation of cell types. 
# This script tries 5 methods, using our scRNA and ATAC from Cusanovich et al. 2018 
# (shendure and furlong lab collab).

wt_scRNA=readRDS( file.path(RESULTS, "object.Rdata"))
results_path = file.path(RESULTS, "../cusanovich_et_al")
dir.create(results_path, showWarnings = F, recursive = T)
cusanovich_et_al_file = file.path(results_path, "cusanovich_et_al.Robj")
merged_file = file.path(results_path, "merged.Robj")

# Load the ATAC data from Furlong / Shendure groups
{
  cusanovich_et_al_counts = read.table(
    "accessory_data/cusanovich_et_al/GSM2706365_10to12.2kbmatrix.txt.gz",
    header = T
  )
  dim(cusanovich_et_al_counts) 
  cusanovich_et_al_features = cusanovich_et_al_counts[ (1:4)]
  cusanovich_et_al_counts   = cusanovich_et_al_counts[-(1:4)]
  cusanovich_et_al_features$i = seq_along(cusanovich_et_al_features$start)
  cusanovich_et_al_cells = data.frame(barcode = colnames(cusanovich_et_al_counts) )
  cusanovich_et_al_cells$j = seq_along(cusanovich_et_al_cells$barcode)
  cusanovich_et_al_celltypes = read.csv("accessory_data/cusanovich_et_al/clusterLabels.csv", header = F) %>% 
    set_colnames(c("cusanovich_cluster", "cusanovich_cell_type"))
  cusanovich_et_al_clusters = read.table("accessory_data/cusanovich_et_al/10to12.densitypeakclusters.txt.gz")
  cusanovich_et_al_tsne = read.table("accessory_data/cusanovich_et_al/10to12.tsnecoords.txt.gz")
  cusanovich_et_al_tsne %<>% set_colnames(c("barcode", "cusanovich_tsne1", "cusanovich_tsne2"))
  cusanovich_et_al_clusters %<>% set_colnames(c("barcode", "cusanovich_cluster"))
  cusanovich_et_al_cells %<>% merge(cusanovich_et_al_tsne)
  cusanovich_et_al_cells %<>% merge(cusanovich_et_al_clusters)
  cusanovich_et_al_cells %<>% merge(cusanovich_et_al_celltypes)
  # These were aligned to dm3, so use that for computing per-gene accessibility scores
  library("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
  drosophila_genes_dm3 = genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
  drosophila_genes_dm3$gene_symbol = convert_fbgn_to_symbol(drosophila_genes_dm3$gene_id)
  drosophila_genes_dm3$gene_name = drosophila_genes_dm3$gene_symbol
  drosophila_genes_dm3$gene_biotype = "gene"
  
  # We don't have a fragments file from GEO or even a tagged BAM. 
  # In lieu of re-processing everything from FASTQ onwards, make a fake fragments file. 
  # Formatted like 10x, with dense-to-sparse heavy lifting by R's Matrix library. 
  # https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments
  fakeFragmentsFile = "accessory_data/cusanovich_et_al/fake_fragments.tsv.gz"
  if(!file.exists(fakeFragmentsFile)){
    cusanovich_et_al_counts %>% 
      as.matrix %>% 
      Matrix::Matrix(sparse = T) %>%
      summary %>%
      merge(cusanovich_et_al_features, by = "i", all.x = T, all.y = F) %>% 
      merge(cusanovich_et_al_cells,    by = "j", all.x = T, all.y = F) %>%
      extract(c("chr", "start", "end", "barcode", "x")) %>%
      set_colnames(c("chrom", "chromStart", "chromEnd", "barcode","readSupport")) %>%
      dplyr::arrange(chrom, chromStart) %>%
      write.table(tools::file_path_sans_ext(fakeFragmentsFile), sep = "\t", col.names = F, row.names = F, quote = F)
    system(paste0("bgzip ", tools::file_path_sans_ext(fakeFragmentsFile)))
    system(paste0("tabix --preset=bed ", fakeFragmentsFile))
  }
}
# Put it in a Seurat object
{
  if(!file.exists(cusanovich_et_al_file)){
    library("Signac")
    atac.assay <- CreateChromatinAssay(
      counts = cusanovich_et_al_counts,
      min.features = 0,
      min.cells = 0,
      max.cells = Inf,
      annotation = drosophila_genes_dm3,
      ranges = GenomicRanges::GRanges(
        seqnames = cusanovich_et_al_features$chr, 
        ranges = IRanges::IRanges(
          cusanovich_et_al_features$start,
          cusanovich_et_al_features$end
        )
      )
    )
    cusanovich_et_al <- CreateSeuratObject(counts = atac.assay, assay = "peaks")
    Fragments(cusanovich_et_al) <- CreateFragmentObject(
      path = fakeFragmentsFile,
      cells = colnames(cusanovich_et_al),
      validate.fragments = FALSE
    )
    
    # Basic unsupervised analysis
    cusanovich_et_al %<>% Signac::FindTopFeatures(min.cutoff = "q90")
    # Exclude X and Y chromosomes to avoid male-female segregation
    VariableFeatures(cusanovich_et_al) %<>% grep(x = ., ignore.case = T, value = T, invert = T, pattern = "chrX|chrY")
    cusanovich_et_al %<>% Signac::RunTFIDF()
    cusanovich_et_al %<>% Signac::RunSVD()
    cusanovich_et_al %<>% Seurat::RunUMAP(reduction = "lsi", dims = 1:30)
    cusanovich_et_al %<>% Seurat::FindNeighbors(reduction = "lsi", dims = 1:30)
    cusanovich_et_al %<>% Seurat::FindClusters(reduction = "lsi", dims = 1:30)
    cusanovich_et_al[["RNA"]] = CreateAssayObject(cusanovich_et_al_gene_activity, min.cells = 0, min.features = 0)
    cusanovich_et_al %<>% NormalizeData(assay = "RNA", normalization_method = "RC")
    {
      pdf(file.path(results_path, "overview.pdf"), width = 7, height = 7)
      DimPlot(cusanovich_et_al, group.by = "seurat_clusters") %>% print
      FeaturePlot(cusanovich_et_al, "nCount_peaks", coord.fixed = T) %>% print
      FeaturePlot(cusanovich_et_al, "nFeature_peaks", coord.fixed = T) %>% print
      salivaryGlandMarkers = 
        read.csv("results/v12/scoring/SalivaryGland/cutoff=2/integratedUMAP/data_driven_markers_best.csv") %>%
        extract2("X") %>%
        intersect(rownames(cusanovich_et_al@assays$RNA))
      for(g in salivaryGlandMarkers){
        try(FeaturePlot(cusanovich_et_al, g, coord.fixed = T) %>% print)
      }
      dev.off()
      saveRDS(cusanovich_et_al, cusanovich_et_al_file)
    } 
  } else {
    cusanovich_et_al = readRDS(cusanovich_et_al_file)
  }
}
cusanovich_et_al_gene_activity = Signac::GeneActivity(cusanovich_et_al, biotypes = NULL)

# RNA to ATAC label transfer with Seurat
if(!"annotationSeuratTransfer" %in% colnames(cusanovich_et_al@meta.data)){
  transfer_anchors = Seurat::FindTransferAnchors(
    reference = wt_scRNA, 
    query = NormalizeData(CreateSeuratObject(cusanovich_et_al_gene_activity))
  )
  cusanovich_et_al[["annotationSeuratTransfer"]] = 
    Seurat::TransferData(transfer_anchors, 
                         refdata = "annotation", 
                         reference = wt_scRNA)$predicted.id
  saveRDS(cusanovich_et_al, cusanovich_et_al_file)
}
# RNA to ATAC label transfer with SingleCellNet
if(!"annotationSCNTransfer" %in% colnames(cusanovich_et_al@meta.data)){
  wt_scRNA@meta.data$barcode = rownames(wt_scRNA@meta.data)
  sharedFeatures = intersect(
    rownames(GetAssayData(wt_scRNA, slot = "scale.data")), 
    rownames(cusanovich_et_al_gene_activity)
  )
  TrainDataSCN = wt_scRNA %>%
    GetAssayData(slot = "scale.data") %>% 
    extract(sharedFeatures, ) 
  SCNClassifier = singleCellNet::scn_train(stTrain  = wt_scRNA@meta.data, 
                                           expTrain = TrainDataSCN,
                                           nTopGenes = 10, 
                                           nRand = 70, 
                                           nTrees = 1000, 
                                           nTopGenePairs = 25,
                                           dLevel = "annotation", 
                                           colName_samp = "barcode")
  classRes_val_all = singleCellNet::scn_predict(
    cnProc=SCNClassifier[['cnProc']], 
    expDat=scale(cusanovich_et_al_gene_activity), 
    nrand = 50
  )
  scn_predictions = classRes_val_all %>% apply(2, function(x) names(which.max(x)))
  cusanovich_et_al[["annotationSCNTransfer"]] = scn_predictions[seq_along(Cells(cusanovich_et_al))]
  saveRDS(cusanovich_et_al, cusanovich_et_al_file)
}

# Process and display the label transfer results
for(labelTransferMethod in c("SCN", "Seurat")){
  cusanovich_et_al %<>% AnnotateClusters(clusterAnnotation, 
                                         by.y = "annotation",
                                         by.x = paste0("annotation", labelTransferMethod, "Transfer"),
                                         namesToAdd = "shortAnnotation",
                                         newNames = paste0("shortAnnotation", labelTransferMethod, "Transfer"))
  DimPlot(cusanovich_et_al, group.by = paste0("shortAnnotation", labelTransferMethod, "Transfer"), label = T) + coord_fixed()
  ggsave(file.path(results_path, paste0("annotation", labelTransferMethod, "Transfer.pdf")), width = 8, height = 5)
}
{
  pdf(file.path(results_path, "comparison.pdf"), width = 5, height = 5)
  X = table(cusanovich_et_al$shortAnnotationSCNTransfer, cusanovich_et_al$shortAnnotationSeuratTransfer) 
  heatmap(X %>% log1p,
          Rowv = NA,
          Colv = NA, 
          xlab = "Seurat",
          ylab = "SCN", 
          scale = "none",
          main = "Log1p total cells with this label")
  dev.off()
}

# Data integration with Seurat
{
  if(file.exists(merged_file)){
    merged_atac_rna = readRDS(merged_file)
  } else {
    anchors = FindIntegrationAnchors(object.list = list(wt_scRNA, CreateSeuratObject(cusanovich_et_al_gene_activity)), 
                                     anchor.features = VariableFeatures(wt_scRNA) )
    merged_atac_rna = IntegrateData(anchorset = anchors)
    merged_atac_rna %<>% ScaleData(assay = "integrated")
    merged_atac_rna %<>% RunPCA( assay = "integrated", reduction.name = "pcaSeurat", npcs = 100)
    merged_atac_rna %<>% RunUMAP( reduction = "pcaSeurat", reduction.name  = "umapSeurat", dims = 1:100 )
    merged_atac_rna %<>% FindNeighbors(reduction = "pcaSeurat")
    merged_atac_rna %<>% FindClusters()
    merged_atac_rna[["clusterSeurat"]] = merged_atac_rna[["seurat_clusters"]]
    
    # This code isn't doing any data integration. I just want to add important metadata to the merged object.
    merged_atac_rna[["cusanovich_cluster"]] = NA
    merged_atac_rna[["cusanovich_cell_type"]] = NA
    merged_atac_rna[["cusanovich_cluster"]][[1]][merged_atac_rna[["modality"]]=="atac"] = 
      cusanovich_et_al[["cusanovich_cluster"]][[1]]
    merged_atac_rna[["cusanovich_cell_type"]][[1]][merged_atac_rna[["modality"]]=="atac"] = 
      cusanovich_et_al[["cusanovich_cell_type"]][[1]]
    merged_atac_rna[["modality"]] = ifelse(is.na(merged_atac_rna[["Phase"]][[1]]), "atac", "rna")
    merged_atac_rna[["barcode_contains"]][[1]] %<>% fillNA(filler = "atac")
    
    # Save for later
    saveRDS(merged_atac_rna, merged_file)
  } 
}

# Data integration with LIGER
if(! "clusterLiger" %in% colnames(merged_atac_rna@meta.data)){
  library(rliger)
  liger <- rliger::createLiger(list(
    rna = wt_scRNA@assays$RNA@counts, 
    atac = cusanovich_et_al_gene_activity
  ))
  liger %<>% rliger::normalize()
  liger %<>% rliger::selectGenes(var.thresh = 0.1)
  liger %<>% rliger::scaleNotCenter()
  liger %<>% rliger::optimizeALS(k=30, use.unshared = TRUE, max_iters =30,thresh=1e-10)
  liger %<>% rliger::quantile_norm()
  liger %<>% rliger::louvainCluster()
  liger %<>% rliger::runUMAP()
  merged_atac_rna[["clusterLiger"]] = liger@clusters
  merged_atac_rna[["umapLiger"]] = SeuratObject::CreateDimReducObject(embeddings = liger@tsne.coords, 
                                                                      key = "umapLiger") 
  saveRDS(merged_atac_rna, merged_file)
}
 
# With Harmony
if(! "clusterHarmony" %in% colnames(merged_atac_rna@meta.data)){
  merged_atac_rna %<>% Seurat::NormalizeData()
  merged_atac_rna %<>% Seurat::FindVariableFeatures()
  Seurat::VariableFeaturePlot(merged_atac_rna)
  merged_atac_rna %<>% Seurat::ScaleData()
  merged_atac_rna %<>% Seurat::RunPCA(npcs = 100)
  merged_atac_rna = harmony::RunHarmony(merged_atac_rna, "barcode_contains", 
                                        plot_convergence = TRUE, 
                                        dims.use = 1:100)
  merged_atac_rna %<>% RunUMAP( reduction = "harmony", reduction.name  = "umapHarmony", dims = 1:100 )
  merged_atac_rna %<>% FindNeighbors( reduction = "harmony", dims = 1:100 )
  merged_atac_rna %<>% FindClusters()
  merged_atac_rna[["clusterHarmony"]] = merged_atac_rna[["seurat_clusters"]]
  saveRDS(merged_atac_rna, merged_file)
}

# Display results
for(integrationMethod in c("Liger", "Seurat", "Harmony")){
  pdf(file.path(results_path, paste0("jointEmbedding", integrationMethod, ".pdf")), width = 8, height = 5)
  DimPlot( merged_atac_rna, group.by = "cusanovich_cell_type", reduction = paste0("umap", integrationMethod), label = T) %>% print
  DimPlot( merged_atac_rna, group.by = "shortAnnotation",      reduction = paste0("umap", integrationMethod), label = T) %>% print
  dev.off()
}



for(feature in c("barcode_contains",
                 "seurat_clusters_atac",
                 "annotations",
                 "modality",
                 "nCount_RNA",
                 "nFeature_RNA",
                 "CG8708",
                 "sage",
                 "PH4alphaSG2")){
  try( DimPlot(    merged_atac_rna, group.by = "clusterSeurat", reduction = "umapSeurat"), silent = T )
  
  try( FeaturePlot(merged_atac_rna, features = feature, coord.fixed = T, reduction = "umapSeurat"), silent = T )
}

