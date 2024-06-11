wt_scRNA=readRDS(file.path("results", ANALYSIS_VERSION, "object.Rdata"))
marker_files = list.files(full.names = T, "accessory_data/v2_2021NOV19")
names(marker_files) = basename(marker_files) %>% gsub(".csv$", "", .)
scores = as.list( names(marker_files) ) %>% setNames(names(marker_files) ) 
i = 1
plots = list()
for(axes in c("umap", "integratedUMAP")){
  for(cell_type in names(marker_files)){
    genes_to_plot = read.csv( marker_files[[cell_type]], header = T) %>% subset(use_in_score = "y")
    genes_to_plot = convert_fbgn_to_symbol(genes_to_plot[["validated_id"]]) %>%
      unique %>%
      intersect(rownames(GetAssayData(wt_scRNA, slot = "counts")))
    cat(cell_type, "genes_to_plot: \n")
    print(genes_to_plot)
    
    # Add up all the genes that we chose as markers
    output_dir_per_celltype = file.path("results", ANALYSIS_VERSION, "scoring", cell_type)
    dir.create(output_dir_per_celltype, recursive = T)
    withr::with_dir(output_dir_per_celltype, {
      scores[[cell_type]] = colSums(GetAssayData(wt_scRNA, slot = "counts")[genes_to_plot, ])
      scores[[cell_type]]  =   scores[[cell_type]]  / wt_scRNA$nCount_RNA 
      scores[[cell_type]] %<>% sqrt
      scores[[cell_type]] %<>% (function(x) (x-mean(x))/sd(x))
      {
        pdf("score_distribution.pdf")
        hist(   scores[[cell_type]], 40, main = paste0("Score for ", cell_type) )
        dev.off()
      }
      for(cutoff in 1:3){
        try({
          current_name = paste0("is_score_above_", cutoff, "_", cell_type)
          wt_scRNA[[current_name]] = ifelse(scores[[cell_type]] > cutoff, cell_type, "other")
          wt_scRNA[[paste0("score_", cell_type)]] = scores[[cell_type]] 
          my_dear = paste0("cutoff=", cutoff, "/", axes)
          dir.create(my_dear, recursive = T)
          withr::with_dir(my_dear, {
            write.table( table( wt_scRNA[[current_name]] ), "num_cells.txt" )
            write.table( table( wt_scRNA[[current_name]][[1]], wt_scRNA$seurat_clusters ), "num_cells.txt" )
            try(silent = T, {
              p = DimPlot(
                wt_scRNA,
                group.by = current_name, 
                cols = c("grey", "blue") %>% setNames(c("other", cell_type)), 
                reduction = axes
              ) + coord_fixed()
              plots[[i]] = p
              i = i + 1
              print(p)
              ggsave("selected_cells.pdf", plot = p, width = 5, height = 5)
            })
            X = Seurat::FindMarkers(wt_scRNA, group.by = current_name, method = "MAST", ident.1 = cell_type)
            X %<>% dplyr::arrange(desc(pct.1/pct.2))
            write.csv(X,           "data_driven_markers.csv" )
            write.csv(head(X, 30), "data_driven_markers_best.csv" )
            Seurat::FindMarkers(wt_scRNA, group.by = current_name, method = "MAST", ident.1 = cell_type,
                                logfc.threshold	= 0, min.pct = 0, min.diff.pct = 0,
                                features = genes_to_plot) %>% 
              dplyr::arrange(p_val) %>% 
              (function(X) X[complete.cases(X), ] ) %>%
              write.csv("handpicked_markers.csv")
          })
        })
      }
    })
  }
}
library(cowplot)
grid = cowplot::plot_grid(plotlist = plots, ncol = 3)
ggsave(paste0("results/", ANALYSIS_VERSION, "/scoring/marker_cutoffs.pdf"), grid, width = 20, height = 40)
ggsave(paste0("results/", ANALYSIS_VERSION, "/scoring/marker_cutoffs.png"), grid, width = 20, height = 40)

grid_strict = cowplot::plot_grid(plotlist = plots[3*(1:(length(plots)/3))], ncol = 4, byrow = T)
ggsave(paste0("results/", ANALYSIS_VERSION, "/scoring/marked_cells_strict_cutoff.pdf"), grid_strict, width = 20, height = 15)
ggsave(paste0("results/", ANALYSIS_VERSION, "/scoring/marked_cells_strict_cutoff.png"), grid_strict, width = 20, height = 15)
saveRDS(wt_scRNA, file.path("results", ANALYSIS_VERSION, "object.Rdata"))
