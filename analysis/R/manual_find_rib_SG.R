library(Seurat)
library(harmony)
library(magrittr)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "manual_find_rib_SG")
dir.create(TARGET_dir)

rib_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_rib", 'manual_celltyping_2_obj.rds'))

rib_object = subset(rib_object, subset = (RNA_snn_res.1 == 0 | RNA_snn_res.1 == 17 | RNA_snn_res.1 == 6))

raw_object = Seurat::CreateSeuratObject(rib_object@assays$RNA@counts, project = 'manual_find_SG')
raw_object@meta.data$batch = rib_object@meta.data$batch
raw_object@meta.data$broad_ct = rib_object@meta.data$manual_celltypes

object = Seurat::NormalizeData(raw_object)
cellCycleMarkers = read.csv("accessory_data/cellCycleMarkers.csv", skip = 1, header = T)
object %<>% CellCycleScoring(s.features = cellCycleMarkers$S.phase.markers., g2m.features = cellCycleMarkers$G2.M.phase.markers.)

genes_in_object = rownames(object@assays$RNA@meta.features)

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

object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
object %<>% ScaleData(features = VariableFeatures(object = object))
 
object@meta.data$broad_ct = rib_object@meta.data$manual_celltypes
#object %<>% RunPCA(features = VariableFeatures(object = object), npcs = 100)
#pc_cutoff = RMTstat::qmp( 1, ndf=length(Cells(object)), pdim=length(VariableFeatures(object)), var=1)
#singular_values = slot(Reductions(object, slot = "pca"), "stdev")
#is_significant = singular_values^2 > pc_cutoff
#num_pc = sum(is_significant)

num_pc = 45 # this is hard coded 
object %<>% RunPCA(features = VariableFeatures(object = object), npcs = num_pc)
object %<>% FindNeighbors(dims = 1:num_pc, k.param = 10) # to really break down the clusters. Probably will need to do the same for the wild types
object %<>% RunUMAP(dim = 1:num_pc)
#object %<>% FindClusters()
pdf(file = file.path(TARGET_dir, "naive_integration.pdf"))
DimPlot(object, group.by = "batch")
dev.off()

pdf(file = file.path(TARGET_dir, "naive_ct.pdf"))
DimPlot(object, group.by = "broad_ct")
dev.off()

set.seed(1)
# harmony integration 
object %<>% harmony::RunHarmony("batch")
object %<>% FindNeighbors(reduction = "harmony", dims = 1:num_pc, k.param = 5)
object %<>% RunUMAP(reduction = "harmony", dim = 1:num_pc)
object %<>% FindClusters(resolution = 2) # super resolution to find the super rare salivary gland cells 

pdf(file = file.path(TARGET_dir, "harmony_integration.pdf"))
DimPlot(object, group.by = "batch")
dev.off()

pdf(file = file.path(TARGET_dir, "harmony_ct.pdf"))
DimPlot(object, group.by = "broad_ct")
dev.off()

FeaturePlot(object, features = 'pip')
object %<>% FindClusters(resolution = 1.5) 

DimPlot(object)
ggsave(filename = file.path(TARGET_dir, 'clustering_epidermis.png'), width = 8, height = 6)

FeaturePlot(object, features = 'pip')
ggsave(filename = file.path(TARGET_dir, 'pip_epidermis.png'), width = 8, height = 6)

FeaturePlot(object, features = 'CG13159')
ggsave(filename = file.path(TARGET_dir, 'CG13159_epidermis.png'), width = 8, height = 6)

VlnPlot(object, features = 'pip', pt.size = 0)
ggsave(filename = file.path(TARGET_dir, 'pip_violin.png'), width = 8, height = 6)

VlnPlot(object, features = 'CG13159', pt.size = 0)
ggsave(filename = file.path(TARGET_dir, 'CG13159_violin.png'), width = 8, height = 6)


DimPlot(object)


object@meta.data$target_cluster = object@meta.data$seurat_clusters == 19 #51
#DimPlot(object)
DimPlot(object,group.by =  'target_cluster')

object@meta.data$finer_annotation = object@meta.data$broad_ct
object@meta.data[object@meta.data$seurat_clusters == 19, 'finer_annotation'] = 'salivary gland'

DimPlot(object, group.by = 'finer_annotation')
ggsave(filename = file.path(TARGET_dir, 'finner_find_SG.png'), width = 8, height = 6)

saveRDS(object, file = file.path(TARGET_dir, 'rib_find_SG_object.rds'))




