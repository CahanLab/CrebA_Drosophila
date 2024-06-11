TARGET_dir = file.path("results", ANALYSIS_VERSION, "Figures/CrebA_cellTypes")
dir.create(TARGET_dir)

##### load in the color palette #####
color_palette = readRDS(file.path('results', ANALYSIS_VERSION, 'ct_color_palettes/ct_color_palette.rds'))

##### load in stage 10-12 CrebA mutants #####
crebA_object = readRDS(file.path('results', ANALYSIS_VERSION, 'manual_annotation_crebA_early/manual_celltype_object.rds'))
DimPlot(crebA_object, group.by = 'manual_celltypes', label = TRUE) + 
  scale_color_manual(values =color_palette) + 
  ggtitle("Stage 10-12 CrebA Mutant Embryos")
ggsave(file.path(TARGET_dir, "stage10-12_embryos_umap.png"), width = 12, height = 8)

DimPlot(crebA_object, group.by = 'manual_celltypes', label = FALSE) + 
  scale_color_manual(values =color_palette) + 
  ggtitle("Stage 10-12 CrebA Mutant Embryos")
ggsave(file.path(TARGET_dir, "stage10-12_embryos_umap_unlabeled.png"), width = 12, height = 8)

##### load in stage 13-16 CrebA mutants #####
crebA_object = readRDS(file.path('results', ANALYSIS_VERSION, 'manual_annotation_crebA/manual_celltype_object.rds'))
DimPlot(crebA_object, group.by = 'manual_celltypes', label = TRUE) + 
  scale_color_manual(values =color_palette) + 
  ggtitle("Stage 13-16 CrebA Mutant Embryos")
ggsave(file.path(TARGET_dir, "stage13-16_embryos_umap.png"), width = 12, height = 8)

DimPlot(crebA_object, group.by = 'manual_celltypes', label = FALSE) + 
  scale_color_manual(values =color_palette) + 
  ggtitle("Stage 13-16 CrebA Mutant Embryos")
ggsave(file.path(TARGET_dir, "stage13-16_embryos_umap_unlabeled.png"), width = 12, height = 8)

##### load in stage 10-12 wildtype #####
wt_object = readRDS(file.path('results', ANALYSIS_VERSION, 'harmonized_wildtype_data/stage10-12_reharmonized_seurat.rds'))
DimPlot(wt_object, group.by = 'new_celltypes', label = TRUE) + 
  scale_color_manual(values =color_palette) + 
  ggtitle("Stage 10-12 Wildtype Embryos")
ggsave(file.path(TARGET_dir, "stage10-12_wt_embryos_umap.png"), width = 12, height = 8)

DimPlot(crebA_object, group.by = 'new_celltypes', label = FALSE) + 
  scale_color_manual(values =color_palette) + 
  ggtitle("Stage 10-12 Wildtype Embryos")
ggsave(file.path(TARGET_dir, "stage10-12_wt_embryos_umap_unlabeled.png"), width = 12, height = 8)

##### load in stage 13-16 widltype #####
wt_object = readRDS(file.path('results', ANALYSIS_VERSION, 'harmonized_wildtype_data/stage13-16_reharmonized_seurat.rds'))
DimPlot(wt_object, group.by = 'new_celltypes', label = TRUE) + 
  scale_color_manual(values =color_palette) + 
  ggtitle("Stage 13-16 WildtypeEmbryos")
ggsave(file.path(TARGET_dir, "stage13-16_wt_embryos_umap.png"), width = 12, height = 8)

DimPlot(wt_object, group.by = 'new_celltypes', label = FALSE) + 
  scale_color_manual(values =color_palette) + 
  ggtitle("Stage 13-16 WildtypeEmbryos")
ggsave(file.path(TARGET_dir, "stage13-16_wt_embryos_umap_unlabeled.png"), width = 12, height = 8)

