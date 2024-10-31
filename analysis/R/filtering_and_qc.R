
# Grab the count data 
filtered_counts = Seurat::Read10X(
  file.path("..", "quantification", CELLRANGER, "filtered_feature_bc_matrix")
)

# Convert flybase to gene symbol
gtf_last_field = unique(read.table("../quantification/reference_genome_info/dmel-all-r6.33.gtf.gz", sep = "\t")[["V9"]])
gene_metadata = data.frame(
  flybase_id =
    gtf_last_field %>%
    strsplit(";") %>%
    sapply(extract2, 1) %>%
    gsub("^.*gene_id ", "", .),
  symbol =
    gtf_last_field %>%
    strsplit(";") %>%
    sapply(extract2, 2) %>%
    gsub("^.*gene_symbol ", "", .)
)
gene_metadata = gene_metadata[!duplicated(gene_metadata),]
table(rownames(filtered_counts)  %in% gene_metadata$flybase_id)
converter = setNames(gene_metadata$symbol, gene_metadata$flybase_id)
rownames(filtered_counts) = converter[rownames(filtered_counts)]

# Make a Seurat object
object = Seurat::CreateSeuratObject(filtered_counts)
object %<>% Seurat::NormalizeData()
# Assemble several quantities useful for QC
# Cell cycle
cellCycleMarkers = read.csv("accessory_data/cellCycleMarkers.csv", skip = 1, header = T)
object %<>% CellCycleScoring(s.features = cellCycleMarkers$S.phase.markers., g2m.features = cellCycleMarkers$G2.M.phase.markers.)

# Certain classes of transcripts
genes_in_object = rownames(object@assays$RNA@meta.features)

# Dan: no need to rerun this. I think there is something wrong with using zcat on MacOS 
#system("zcat ../../quantification/reference_genome_info/dmel-all-r6.33.gtf.gz | grep mitochondrion | cut -f9 | cut -d';' -f1 | cut -f2 -d' ' | sort | uniq > accessory_data/mitochondrially_encoded_genes.tsv")
#system("zcat ../../quantification/reference_genome_info/dmel-all-r6.33.gtf.gz | grep rDNA | cut -f9 | cut -d';' -f1 | cut -f2 -d' ' | sort | uniq > accessory_data/ribosomal_rna_genes.tsv")

mitochondrially_encoded_genes = read.table("accessory_data/mitochondrially_encoded_genes.tsv", header = F)[[1]]
ribosomal_rna_genes           = read.table("accessory_data/ribosomal_rna_genes.tsv", header = F)[[1]]
ribosomal_protein_genes       = grep("^rpl|^rps", ignore.case = T, value = T, genes_in_object)
object[["mitochondrial_transcript_total_expression"]] =
  GetAssayData(object, "counts") %>%
  extract( convert_fbgn_to_symbol( mitochondrially_encoded_genes ), ) %>%
  colSums
object[["ribosomal_rna_total_expression"]] =
  GetAssayData(object, "counts") %>%
  extract( convert_fbgn_to_symbol( ribosomal_rna_genes ), ) %>%
  colSums
object[["ribosomal_protein_total_expression"]] =
  GetAssayData(object, "counts") %>%
  extract( ribosomal_protein_genes, ) %>%
  colSums
object[["mitochondrial_transcript_norm_expression"]] = object[["mitochondrial_transcript_total_expression"]] / object$nCount_RNA
object[["ribosomal_rna_norm_expression"           ]] = object[["ribosomal_rna_total_expression"           ]] / object$nCount_RNA
object[["ribosomal_protein_norm_expression"       ]] = object[["ribosomal_protein_total_expression"       ]] / object$nCount_RNA
object[["log10_nCount_RNA"]] = object[["nCount_RNA"]] %>% log10


# Basic analysis, determining number of PC's with Marchenko-Pastur upper bound
# Dan: I don't think you need to normalize the data again 
# object %<>% Seurat::NormalizeData() 

object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
object %<>% ScaleData(features = VariableFeatures(object = object))
object %<>% RunPCA(features = VariableFeatures(object = object), npcs = 100)

pc_cutoff = RMTstat::qmp( 1, ndf=length(Cells(object)), pdim=length(VariableFeatures(object)), var=1)
singular_values = slot(Reductions(object, slot = "pca"), "stdev")
is_significant = singular_values^2 > pc_cutoff
num_pc = sum(is_significant)
object %<>% FindNeighbors()
object %<>% FindClusters()
object %<>% RunUMAP(dims = 1:num_pc)

# Show results before removing any cells, then start removing cells
dir.create(file.path(RESULTS, "before_cleaning"), recursive = T, showWarnings = F)
withr::with_dir(file.path(RESULTS, "before_cleaning"), {
  write.csv( cbind( object[["seurat_clusters"]], object[["umap"]]@cell.embeddings ), 
             "umap_and_clustering.csv")
  pdf("pairwise_qc.pdf", width = 12, height = 10)
  pairs(object@meta.data[c("mitochondrial_transcript_norm_expression", 
                           "ribosomal_protein_norm_expression", 
                           "log10_nCount_RNA")], 
        col = COLORSCALES$phase[object@meta.data$Phase],
        pch = ".") 
  par(xpd=T)
  legend("topleft", legend = names(COLORSCALES$phase), fill = COLORSCALES$phase)
  dev.off()
  MakeBasicPlots(object, reduction = "umap")
  
  # Remove any barcodes with high mitochondrial transcript counts 
  write.csv( table(removed = object$mitochondrial_transcript_norm_expression>=0.25),"mitochondria_removal.csv")
  object %<>% subset( mitochondrial_transcript_norm_expression <= 0.25)
  
  # Remove any cluster with few marker genes 
  # on suspicion that it contains mostly empties.
  if('markers.csv' %in% list.files()) { 
    markers = read.csv("markers.csv", row.names = 1)
  } else { 
    markers = FindAllMarkers(object)
    write.csv(markers, "markers.csv")
  }

  how_many_markers = 
    markers %>% 
    subset(avg_log2FC > 0.25) %>%
    extract2("cluster") %>%
    table(cluster = .) %>% 
    as.data.frame() %>% 
    set_colnames(c("cluster", "marker count")) 
  clusters_remove = 
    subset(how_many_markers, `marker count` <= 50) %>%
    extract2("cluster")
  table(cluster = object$seurat_clusters) %>%
    as.data.frame() %>%
    set_colnames(c("cluster", "cell count")) %>%
    dplyr::mutate(is_removed = cluster %in% clusters_remove) %>%
    merge(how_many_markers) %>%
    write.csv( "cluster_stats.csv" )
  object = object[,!(object$seurat_clusters %in% clusters_remove)]
  # Filter out doublets after removing empties
  # 10X doublet rate is 0.9% per 1,000 cells
  n_cell = length(Cells(object))
  dub_proportion = 0.009*n_cell/1000
  dubbitydub = DoubletFinder::doubletFinder_v3(object, 
                                               PCs = 1:100, 
                                               nExp = round(dub_proportion*n_cell),
                                               pK = 100 / n_cell)
  stopifnot( length( unique( dubbitydub@meta.data$pANN ) ) > 10 )
  object$DoubletFinder_score = dubbitydub@meta.data$pANN
  object$DoubletFinder_call = dubbitydub@meta.data$DF.classifications
  rm(dubbitydub); gc()
  write.csv( table(object$DoubletFinder_call), "doublet_removal.csv")
  object %<>% subset(DoubletFinder_call == "Singlet")
})

# Save analysis results
withr::with_dir(RESULTS, {
  saveRDS(object, "object.Rdata")
})



