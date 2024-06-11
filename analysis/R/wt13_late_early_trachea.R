
TARGET_dir = file.path("results", ANALYSIS_VERSION, "wt_late_early_trachea")
dir.create(TARGET_dir)

wt_late_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object3.rds"))

##########################
# TODO the below will be changed as time goes on 
wt_early_object = readRDS(file.path("results", ANALYSIS_VERSION, "early_wt12_integrated/BDGP_automated_annotation_object.rds"))
sub_wt_early_object = subset(wt_early_object, subset = Integrated_tentativeCellType == 'tracheal primordium')
sub_wt_early_object$experimental_condition = 'early'
########################

sub_wt_late_object = subset(wt_late_object, subset = manual_celltypes == "Trachea")
sub_wt_late_object$experimental_condition = 'late'

combined_ct_object = merge(sub_wt_early_object, sub_wt_late_object)
object = Seurat::CreateSeuratObject(combined_ct_object@assays$RNA@counts, project = 'Trachea')
object@meta.data$experimental_condition = combined_ct_object@meta.data$experimental_condition
object@meta.data$batch = paste0(combined_ct_object@meta.data$experimental_condition, "_", combined_ct_object@meta.data$batch)

# now start analysis 
object = Seurat::NormalizeData(object)

object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
object %<>% ScaleData(features = VariableFeatures(object = object))
object %<>% RunPCA(features = VariableFeatures(object = object), npcs = 100)

pc_cutoff = RMTstat::qmp( 1, ndf=length(Cells(object)), pdim=length(VariableFeatures(object)), var=1)
singular_values = slot(Reductions(object, slot = "pca"), "stdev")
is_significant = singular_values^2 > pc_cutoff
num_pc = sum(is_significant)

object %<>% RunPCA(features = VariableFeatures(object = object), npcs = num_pc)
object %<>% FindNeighbors(dims = 1:num_pc)
object %<>% RunUMAP(dim = 1:num_pc)
object %<>% FindClusters()

DimPlot(object, group.by = 'batch')

DimPlot(object)
ggsave(filename = file.path(TARGET_dir, 'naive_cluster.png'), width = 8, height = 6)

DimPlot(object, group.by = 'experimental_condition')
ggsave(filename = file.path(TARGET_dir, 'naive_experimental_condition.png'), width = 8, height = 6)

DimPlot(object, group.by = 'batch')
ggsave(filename = file.path(TARGET_dir, 'naive_batch.png'), width = 8, height = 6)
saveRDS(object, file = file.path(TARGET_dir, "ct_object.rds"))

# this is with harmony integration 
object %<>% harmony::RunHarmony("batch")
object %<>% FindNeighbors(reduction = "harmony", dims = 1:num_pc)
object %<>% RunUMAP(reduction = "harmony", dim = 1:num_pc)
object %<>% FindClusters()

DimPlot(object)
ggsave(filename = file.path(TARGET_dir, 'harmony_cluster.png'), width = 8, height = 6)

DimPlot(object, group.by = 'experimental_condition')
ggsave(filename = file.path(TARGET_dir, 'harmony_experimental_condition.png'), width = 8, height = 6)

DimPlot(object, group.by = 'batch')
ggsave(filename = file.path(TARGET_dir, 'harmony_batch.png'), width = 8, height = 6)
saveRDS(object, file = file.path(TARGET_dir, "ct_object_harmony.rds"))


########################
library(monocle3)
seurat_object = readRDS(file.path(TARGET_dir, "ct_object.rds"))
expression_matrix = seurat_object@assays$RNA@counts
cell_metadata = seurat_object@meta.data
gene_annotation = seurat_object@assays$RNA@meta.features
gene_annotation$gene_short_name = rownames(gene_annotation)
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)
cds <- align_cds(cds, alignment_group = "batch")

cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "batch", cell_size = 1, label_cell_groups = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_batch_corrected_UMAP.png'), width = 8, height = 6)

cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE, cell_size = 1)
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "batch",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, cell_size = 1)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1, show_trajectory_graph = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_batch_corrected_pt.png'), width = 8, height = 6)
saveRDS(cds, file = file.path(TARGET_dir, "monocle3_batch_correct_object.rds"))

cds = readRDS(file.path(TARGET_dir, "monocle3_batch_correct_object.rds"))
plot_cells(cds,
           gene = c("Osi17"),
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1, show_trajectory_graph = FALSE)

##############
seurat_object = readRDS(file.path(TARGET_dir, "ct_object.rds"))

# you can uncomment the below comment to get the full dataset 
# seurat_object = subset(seurat_object, subset = batch != 'early_rep1')

expression_matrix = seurat_object@assays$RNA@counts
cell_metadata = seurat_object@meta.data
gene_annotation = seurat_object@assays$RNA@meta.features
gene_annotation$gene_short_name = rownames(gene_annotation)
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)

cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "batch", cell_size = 1, label_cell_groups = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_UMAP.png'), width = 8, height = 6)

plot_cells(cds,
           genes=c("toe"),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE, cell_size = 2)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_toe.png'), width = 8, height = 6)

cds <- cluster_cells(cds, resolution = 1e-2)
plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE, cell_size = 1)
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "batch",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1, show_trajectory_graph = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_pt.png'), width = 8, height = 6)
saveRDS(cds, file = file.path(TARGET_dir, "monocle3_no_batch_correct_object.rds"))
