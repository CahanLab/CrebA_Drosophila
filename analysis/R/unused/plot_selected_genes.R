source("R/set_up_environment.R")
wt_scRNA=readRDS("results/v3/initial_exploration.Rdata")
marker_files = list.files(full.names = T, "accessory_data/v2_2021NOV19")
names(marker_files) = basename(marker_files)  %>% gsub("\\.csv$", "", .)
for(cell_type in names(marker_files)){
  genes_to_plot = read.csv( marker_files[[cell_type]], header = T)
  genes_to_plot = convert_fbgn_to_symbol(genes_to_plot[["validated_id"]])
  cat(cell_type, "genes_to_plot: \n")
  print(genes_to_plot)
  # Show UMAP, clusters, per-cell depth, and manually selected genes
  {
    pdf(paste0("results/v3/", cell_type, ".pdf"), width = 10, height = 8)
    DimPlot(wt_scRNA, group.by = "seurat_clusters", label = T) %>% print
    for( gene in genes_to_plot ){
      try({FeaturePlot(wt_scRNA, features = gene, coord.fixed = T) %>% print})
    }
    dev.off()
  }
}
