object=readRDS(file.path(RESULTS, "object.Rdata"))
Idents(object)=object[["seurat_clusters"]][[1]]

# Crank out cluster markers
Idents(object)=object[["seurat_clusters"]][[1]]
differentially_expressed = FindAllMarkers(object, test.use = "bimod")
withr::with_dir(
  file.path(RESULTS), 
  {
    write.csv(differentially_expressed, "marker_genes.csv")
    differentially_expressed %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n=50, wt = avg_log2FC) %>%
      dplyr::arrange(cluster) %>%
      write.csv("top_marker_genes_1-40.csv")
  }
)

# Functional enrichment analysis
library("enrichR")
enrichR::setEnrichrSite("FlyEnrichr")
enrichR::listEnrichrDbs() 
doEnrichr = function(cluster_markers){
  cluster_markers %>% 
    extract2("gene") %>%
    enrichR::enrichr(
      genes = ., 
      databases = c(
        "GO_Molecular_Function_2018",
        "GO_Biological_Process_2018"
      )
    ) %>% 
    dplyr::bind_rows(.id = "database") %>%
    dplyr::group_by(database) %>% 
    dplyr::top_n(wt = -Adjusted.P.value, n = 5) 
}
enrichment_results = 
  differentially_expressed %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n=50, wt = avg_log2FC) %>%
  as(., "data.table") %>%
  split(by="cluster") %>% 
  lapply(doEnrichr) %>% 
  dplyr::bind_rows(.id = "cluster") 
withr::with_dir(
  file.path(RESULTS), 
  {
    write.csv(enrichment_results, "enrichment_results.csv")
  }
)

# Use previously annotated marker genes for a rough guess at cell type identification.
annotationPerGene    = paste0("accessory_data/annotations/2022MAY13/genes/rib_rep1.csv") # rhis part was modified to include rib_rep1
if(file.exists(annotationPerGene)){
  geneAnnotation = read.csv(annotationPerGene)
  geneAnnotation %<>% extract(c("gene", "annotation"))
  # Guess cluster: max overlap from prior markers
  clusterAnnotation = 
    differentially_expressed %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n=20, wt = avg_log2FC) %>%
    merge(., geneAnnotation, all.x = T, all.y = F, by = "gene") %>%
    dplyr::arrange(cluster, avg_log2FC) %>% 
    with(table(cluster, annotation)) %>% 
    as.data.frame %>% 
    dplyr::rename(n_supporting_genes = Freq) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(n=1, order_by = n_supporting_genes, with_ties = T) %>%
    dplyr::arrange(cluster)
  clusterAnnotation$annotation %<>% as.character()
  # Deal with low overlap
  clusterAnnotation$annotation[clusterAnnotation$n_supporting_genes<=2] = "unknown"
  # Deal with ties
  clusterAnnotation %<>% 
    group_by(cluster) %>%
    dplyr::summarise(n_supporting_genes = mean(n_supporting_genes), 
                     annotation = paste0(unique(annotation), collapse = " or\n"))
  # Record results
  write.csv(clusterAnnotation, file.path(RESULTS, "tentativeCellTypes.csv"))
  object %<>% AnnotateClusters(clusterAnnotation[c("cluster", "annotation")], 
                                 namesToAdd = "annotation", 
                                 newNames = "tentativeCellType", 
                                   by.x = "seurat_clusters", by.y = "cluster")
  DimPlot(object, group.by = "tentativeCellType", label = T, label.size = 5) +
    coord_fixed() + 
    ggtitle("Tentative cell type")
  ggsave( file.path( RESULTS, "tentativeCellTypes.pdf" ), width = 10, height = 7 )
}

