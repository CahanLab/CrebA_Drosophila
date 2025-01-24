TARGET_dir = file.path("results", ANALYSIS_VERSION, "spearman_sim_early_creba_wt")
dir.create(TARGET_dir)

##### load in the data sets ######
wt_object = readRDS(file.path('results', ANALYSIS_VERSION, "harmonized_wildtype_data/stage10-12_reharmonized_seurat.rds"))
wt_object@meta.data$manual_celltypes = wt_object@meta.data$new_celltypes

mut_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_crebA_early/manual_celltype_object.rds"))

calc_mean_spearman <- function(sub_wt_obj, sub_mut_obj) {
  igenes = intersect(rownames(sub_wt_obj@assays$RNA@data), rownames(sub_mut_obj@assays$RNA@data))
  sub_wt_exp = sub_wt_obj@assays$RNA@data[igenes, ]
  sub_mut_exp = sub_mut_obj@assays$RNA@data[igenes, ]
  spearman_corr = cor(apply(sub_wt_exp, MARGIN = 1, FUN = mean), apply(sub_mut_exp, MARGIN = 1, FUN = mean))
  return(spearman_corr)
}

intersect_cts = intersect(wt_object@meta.data$manual_celltypes, mut_object@meta.data$manual_celltypes) 
plot_df = data.frame(celltypes = intersect_cts)
plot_df$creba_exp = NA
plot_df$spearman_corr = NA
rownames(plot_df) = plot_df$celltypes

for(tmp_ct in intersect_cts) {
  sub_wt_obj = subset(wt_object, subset = manual_celltypes == tmp_ct) 
  sub_mut_obj = subset(mut_object, subset = manual_celltypes == tmp_ct)
  
  creba_exp = mean(sub_mut_obj@assays$RNA@data['CrebA', ])
  spearman_corr = calc_mean_spearman(sub_wt_obj, sub_mut_obj)
  plot_df[tmp_ct, 'creba_exp'] = creba_exp
  plot_df[tmp_ct, 'spearman_corr'] = spearman_corr
}

p = ggplot(plot_df, aes(x = creba_exp, y = spearman_corr)) +
  geom_point(color = "blue", size = 3) +  # Plot points
  geom_text(aes(label = celltypes), vjust = -0.5, hjust = 0.5) +  # Label points
  labs(title = "Spearman correlation similarity vs CrebA exp",
       x = "CrebA expression",
       y = "Spearman correlation between mutant and wildtype") +
  xlim(c(-0.2, 1.7)) +
  theme_cowplot()

ggsave(filename = file.path(TARGET_dir, 'spearman_corr.png'), plot = p, width = 12, height = 8)



