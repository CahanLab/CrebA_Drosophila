TARGET_dir = file.path("results", ANALYSIS_VERSION, "check_depth")
dir.create(TARGET_dir)

wt_object = readRDS(file.path('results', ANALYSIS_VERSION, "harmonized_wildtype_data/stage10-12_reharmonized_seurat.rds"))
mut_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_crebA_early/manual_celltype_object.rds"))

# generate the count dataframe 
wt_counts_df = data.frame('n_counts' = wt_object@meta.data$nCount_RNA, 
                          'n_features' = wt_object@meta.data$nFeature_RNA, 
                          'log10_n_counts' = wt_object@meta.data$log10_nCount_RNA,
                          'batch' = paste0('wt-', wt_object@meta.data$batch), 
                          'condition' = 'wt', 
                          'celltype' = wt_object@meta.data$new_celltypes)

mut_counts_df = data.frame('n_counts' = mut_object@meta.data$nCount_RNA, 
                          'n_features' = mut_object@meta.data$nFeature_RNA, 
                          'log10_n_counts' = mut_object@meta.data$log10_nCount_RNA,
                          'batch' = paste0('mut-', mut_object@meta.data$batch), 
                          'condition' = 'mut', 
                          'celltype' = mut_object@meta.data$manual_celltypes)

combined_counts_df = rbind(wt_counts_df, mut_counts_df)

p = ggplot(combined_counts_df, aes(x = batch, y = log10_n_counts)) +
  geom_violin(fill = "skyblue", alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +  # optional: overlay boxplot
  theme_minimal() +
  labs(title = "Read depth of early wild-type and CrebA mutant samples",
       x = "Batch",
       y = "log10(n_counts)") + 
  theme(plot.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14))
ggsave(file.path(TARGET_dir, 'read_depth.png'), plot = p, width = 8, height = 5)

# look at the 
p = ggplot(combined_counts_df, aes(x = batch, y = log10_n_counts)) +
  geom_violin(fill = "skyblue", alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA) +  # optional: overlay boxplot
  theme_minimal() +
  facet_wrap(~ celltype) +  # facet by condition
  labs(title = "Read depth of early wild-type and CrebA mutant samples",
       x = "Batch",
       y = "log10(n_counts)") + 
  theme(plot.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        strip.text = element_text(size = 12))
ggsave(file.path(TARGET_dir, 'read_depth_celltypes.png'), plot = p, width = 16, height = 12)

