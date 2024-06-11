# Compare cell type diversity between cells, nuclei, and unknowns. Do they have 
# the same set of stuff, just damaged or stripped down in various ways?
# Or do they have different stuff?
withr::with_dir(file.path("results", ANALYSIS_VERSION), {
  objects = list(
    cells = readRDS("cells/object.Rdata"),
    nuclei = readRDS("nuclei/object.Rdata"),
    unknown = readRDS("neither/object.Rdata")
  )
})

classifiers = lapply(objects, BuildClassifier)
withr::with_dir(file.path("results", ANALYSIS_VERSION), {
  saveRDS(classifiers, "classifiers.Rdata")
})
for(biotype_reference in names(objects)){
  for(biotype_query in names(objects)){
    objects[[biotype_query]][[paste0("label_from_",biotype_reference)]] = 
      ApplyClassifier(objects[[biotype_query]], classifiers[[biotype_reference]])
  }
}

frequencies = list(); i = 0
for(biotype_reference in names(objects)){
  for(biotype_query in names(objects)){
    i = i + 1
    frequencies[[i]] = data.frame(table(assigned_cluster = objects[[biotype_query]][[paste0("label_from_",biotype_reference)]] ))
    frequencies[[i]][["Freq"]] %<>% (function(x) round((x/sum(x))*100, digits = 1))
    frequencies[[i]][["biotype_query"]] = biotype_query
    frequencies[[i]][["biotype_reference"]] = biotype_reference
  }
}
frequencies %<>% data.table::rbindlist()
frequencies %<>% tidyr::pivot_wider(names_from = "biotype_query",
                                    id_cols = c("assigned_cluster", "biotype_reference"),
                                    values_from = c("Freq"))
withr::with_dir(file.path("results", ANALYSIS_VERSION), {
  # Cells vs nuclei
  ggplot(subset(frequencies, biotype_reference == "cells")) + 
    geom_text(aes(cells, nuclei, label = assigned_cluster))  + 
    geom_abline(slope = 1) + 
    ggtitle("Cluster abundance when classifying nuclei based on cells")
  ggsave("classify_nuclei_based_on_cells.pdf", height = 5, width = 5)
  ggplot(subset(frequencies, biotype_reference == "nuclei")) + 
    geom_text(aes(cells, nuclei, label = assigned_cluster))  + 
    geom_abline(slope = 1) + 
    ggtitle("Cluster abundance when classifying cells based on nuclei")
  ggsave("classify_cells_based_on_nuclei.pdf", height = 5, width = 5)
  DimPlot(objects$cells, group.by = "label_from_nuclei", label = T) + 
    ggtitle("Cells classified into nuclei-defined clusters")
  ggsave("Cells classified into nuclei-defined clusters.pdf", width = 5, height = 5)
  DimPlot(objects$nuclei, group.by = "label_from_cells", label = T) + 
    ggtitle("Nuclei classified into cell-defined clusters")
  ggsave("Nuclei classified into cell-defined clusters.pdf", width = 5, height = 5)
  
  # Cells vs unknown
  ggplot(subset(frequencies, biotype_reference == "cells")) + 
    geom_text(aes(cells, unknown, label = assigned_cluster))  + 
    geom_abline(slope = 1) + 
    ggtitle("Cluster abundance when classifying unknowns based on cells")
  ggsave("classify_unknown_based_on_cells.pdf", height = 5, width = 5)
  ggplot(subset(frequencies, biotype_reference == "unknown")) + 
    geom_text(aes(cells, unknown, label = assigned_cluster))  + 
    geom_abline(slope = 1) + 
    ggtitle("Cluster abundance when classifying cells based on unknowns")
  ggsave("classify_cells_based_on_unknown.pdf", height = 5, width = 5)
  
  DimPlot(objects$cells, group.by = "label_from_unknown", label = T) + 
    ggtitle("Cells classified into unknown-defined clusters")
  ggsave("Cells classified into unknown-defined clusters.pdf", width = 5, height = 5)
  DimPlot(objects$unknown, group.by = "label_from_cells", label = T) + 
    ggtitle("Unknowns classified into cell-defined clusters")
  ggsave("Unknowns classified into cell-defined clusters.pdf", width = 5, height = 5)
})
