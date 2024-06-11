
TARGET_dir = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_germ")
dir.create(TARGET_dir)

wt_late_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object3.rds"))

##########################
# TODO the below will be changed as time goes on 
wt_early_object = readRDS(file.path("results", ANALYSIS_VERSION, "early_wt12_integrated/BDGP_automated_annotation_object.rds"))
sub_wt_early_object = subset(wt_early_object, subset = Integrated_tentativeCellType == 'germ cell')
sub_wt_early_object$experimental_condition = 'early'
########################

sub_wt_late_object = subset(wt_late_object, subset = manual_celltypes == "Germ Cell")
sub_wt_late_object$experimental_condition = 'late'

combined_ct_object = merge(sub_wt_early_object, sub_wt_late_object)
object = Seurat::CreateSeuratObject(combined_ct_object@assays$RNA@counts, project = 'Germ Cell')
object@meta.data$experimental_condition = combined_ct_object@meta.data$experimental_condition
object@meta.data$batch = paste0(combined_ct_object@meta.data$experimental_condition, "_", combined_ct_object@meta.data$batch)

object = Seurat::NormalizeData(object)
cellCycleMarkers = read.csv("accessory_data/cellCycleMarkers.csv", skip = 1, header = T)
object %<>% CellCycleScoring(s.features = cellCycleMarkers$S.phase.markers., g2m.features = cellCycleMarkers$G2.M.phase.markers.)

# start monocle3 
expression_matrix = object@assays$RNA@counts
cell_metadata = object@meta.data
gene_annotation = object@assays$RNA@meta.features
gene_annotation$gene_short_name = rownames(gene_annotation)
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)

cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "batch", cell_size = 1, label_cell_groups = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_UMAP.png'), width = 8, height = 6)

cds <- cluster_cells(cds, resolution = 1e-4)
plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE, cell_size = 1)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_cluster.png'), width = 8, height = 6)

marker_test_res <- monocle3::top_markers(cds, group_cells_by="cluster", 
                                         reference_cells=1000, cores=8)

write.csv(marker_test_res, file = file.path(TARGET_dir, "top_markers_monocle.csv"))
rank_sum_results = presto::wilcoxauc(normalized_counts(cds), cds@clusters$UMAP$clusters)
write.csv(rank_sum_results, file = file.path(TARGET_dir, "rank_sum_test.csv"))

########################################
# let's solve the identity of cluster 1 first 
dir.create(file.path(TARGET_dir, 'cluster_process'))
marker_test_res = read.csv(file.path(TARGET_dir, "top_markers_monocle.csv"), row.names = 1)
withr::with_dir(file.path(TARGET_dir, 'cluster_process'), {
  for(cluster in unique(marker_test_res$cell_group)) {
    sub_marker_test_res = marker_test_res[marker_test_res$cell_group == cluster, ]
    enrichment_results = enrichR::enrichr(
      genes = sub_marker_test_res$gene_id, 
      databases = c(
        "GO_Biological_Process_2018", 
        "GO_Molecular_Function_2018"
      )
    )
    biological_analysis = enrichment_results$GO_Biological_Process_2018
    write.csv(biological_analysis, file = paste0(cluster, "_biological_process.csv"))
  }
})

plot_cells(cds,
           genes=c("nur", 'RpL38', 'RpL30'),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE, cell_size = 2)
plot_cells(cds, color_cells_by = "Phase", label_cell_groups = FALSE, cell_size = 1, show_trajectory_graph = FALSE)

plot_cells(cds, color_cells_by = "nFeature_RNA", label_cell_groups = FALSE, cell_size = 1, show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by = "nCount_RNA", label_cell_groups = FALSE, cell_size = 1, show_trajectory_graph = FALSE)

cds@colData$log_nCount = log(cds@colData$nCount_RNA)
plot_cells(cds, color_cells_by = "log_nCount", label_cell_groups = FALSE, cell_size = 1, show_trajectory_graph = FALSE)

cds@colData$log_nFeature = log(cds@colData$nFeature_RNA)
plot_cells(cds, color_cells_by = "log_nFeature", label_cell_groups = FALSE, cell_size = 1, show_trajectory_graph = FALSE)

ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_phase.png'), width = 8, height = 6)

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
