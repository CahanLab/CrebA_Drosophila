TARGET_dir = file.path("results", ANALYSIS_VERSION, "Figures/CrebA_stats")
dir.create(TARGET_dir)

##### load in stage 10-12 CrebA mutants #####
color_palette = readRDS(file.path('results', ANALYSIS_VERSION, 'ct_color_palettes/ct_color_palette.rds'))
crebA_object = readRDS(file.path('results', ANALYSIS_VERSION, 'manual_annotation_crebA_early/manual_celltype_object.rds'))

# make the batch UMAPs 
DimPlot(crebA_object, group.by = 'batch') + 
  ggtitle("Stage 10-12 CrebA Mutant Embryos")
ggsave(file.path(TARGET_dir, "batch_creba_umap.png"), width = 10, height = 7)

# make the seurat clusters 
DimPlot(crebA_object, group.by = 'seurat_clusters') + 
  ggtitle("Stage 10-12 CrebA Mutant Embryos")
ggsave(file.path(TARGET_dir, "clusters_creba_umap.png"), width = 10, height = 7)

# make the sample table 
samp_tab = crebA_object@meta.data

samp_tab %>% 
  group_by(batch) %>% 
  summarise(mean_nCount = mean(nCount_RNA))

samp_tab %>% 
  group_by(batch) %>% 
  summarise(median_nCount = median(nCount_RNA))

samp_tab %>% 
  group_by(batch) %>% 
  summarise(mean_nCount = mean(nFeature_RNA))

samp_tab %>% 
  group_by(batch) %>% 
  summarise(median_nCount = median(nFeature_RNA))

# make the violin plot 
VlnPlot(crebA_object, features = 'nCount_RNA', group.by = 'batch', pt.size = 0)
