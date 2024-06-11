object=readRDS(file.path(RESULTS, "object.Rdata"))
Idents(object)=object[["integratedClusters"]][[1]]

# Retrieve the latest manual annotations
annotationPerCluster = paste0("accessory_data/annotations/2022MAY10/", SAMPLE, ".csv")
clusterAnnotation = read.csv(annotationPerCluster)
clusterAnnotation$annotation %<>% tolower

# Add them to the object and plot + save
withr::with_dir(
  file.path(RESULTS), 
  {
    object %<>% AnnotateClusters(clusterAnnotation, by.x = "integratedClusters", by.y = "cluster")
    DimPlot(object, group.by = "shortAnnotation", label = T, label.size = 5) +
      coord_fixed() + 
      scale_color_manual(
        labels = setNames(
          paste0(clusterAnnotation$annotation, " (", clusterAnnotation$shortAnnotation, ")"),
          clusterAnnotation$shortAnnotation
        ),
        values = scales::hue_pal()(length(clusterAnnotation$annotation))
      ) +
      ggtitle("Cell type")
    ggsave("cellTypes.pdf", width = 14, height = 7)
    save3DUMAPCategorical(object, umapkey = "umap3D", group.by = "shortAnnotation")
    table(object[["annotation"]]) %>%
      as.data.frame() %>% 
      set_colnames(c("CellType", "Count")) %>%
      dplyr::mutate(Percent = round(digits = 1, 100*Count/sum(Count))) %>%
      write.table("cellTypeCounts.csv")
    saveRDS(object, "object.Rdata")
    # Also save celltype annotations alone in a form that is small enough to check into version control. 
    write.csv(object[["annotation"]], "annotations_by_barcode.csv") 
  }
)


