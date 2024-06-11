library(ggplot2)
library(Seurat)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "Chitin_Comparison")
object = readRDS(file.path("results", ANALYSIS_VERSION, "wt_123_compare/all_3/harmony_object.rds"))

DimPlot(object, group.by = 'batch')
ggsave(file.path("results", ANALYSIS_VERSION, "Chitin_Comparison", paste0("batch_plot.png")), width = 8, height = 6)

gene = 'Cpr65Ea'
FeaturePlot(object, features = gene)
ggsave(file.path("results", ANALYSIS_VERSION, "Chitin_Comparison", paste0(gene, "_feature.png")), width = 8, height = 6)

gene = 'Cpr65Eb'
FeaturePlot(object, features = gene)
ggsave(file.path("results", ANALYSIS_VERSION, "Chitin_Comparison", paste0(gene, "_feature.png")), width = 8, height = 6)

gene = 'Cpr64Ab'
FeaturePlot(object, features = gene)
ggsave(file.path("results", ANALYSIS_VERSION, "Chitin_Comparison", paste0(gene, "_feature.png")), width = 8, height = 6)

gene = 'Ccp84Ae'
FeaturePlot(object, features = gene)
ggsave(file.path("results", ANALYSIS_VERSION, "Chitin_Comparison", paste0(gene, "_feature.png")), width = 8, height = 6)

gene = 'Cpr64Ad'
FeaturePlot(object, features = gene)
ggsave(file.path("results", ANALYSIS_VERSION, "Chitin_Comparison", paste0(gene, "_feature.png")), width = 8, height = 6)

gene = 'Ccp84Ag'
FeaturePlot(object, features = gene)
ggsave(file.path("results", ANALYSIS_VERSION, "Chitin_Comparison", paste0(gene, "_feature.png")), width = 8, height = 6)
