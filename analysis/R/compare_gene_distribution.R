TARGET_dir = file.path("results", ANALYSIS_VERSION, "compare_gene_distribution")
dir.create(TARGET_dir)

wt_object = readRDS(file.path('results', ANALYSIS_VERSION, "harmonized_wildtype_data/stage10-12_reharmonized_seurat.rds"))
wt_object@meta.data$manual_celltypes = wt_object@meta.data$new_celltypes

mut_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_crebA_early/manual_celltype_object.rds"))

common_celltypes = intersect(wt_object@meta.data$manual_celltypes, mut_object@meta.data$manual_celltypes)

##### look at salivary gland cells #####
celltype = 'Salivary Gland'

sub_wt_obj = subset(wt_object, subset = manual_celltypes == celltype) 
sub_wt_obj@meta.data$experimental_condition = 'Wt'

sub_mut_obj = subset(mut_object, subset = manual_celltypes == celltype) 
sub_mut_obj@meta.data$experimental_condition = 'mut'

combined_obj = merge(sub_wt_obj, sub_mut_obj)
object = Seurat::CreateSeuratObject(combined_obj@assays$RNA@counts, project = 'find_diff')
object@meta.data$experimental_condition = combined_obj@meta.data$experimental_condition
object@meta.data$batch = paste0(combined_obj@meta.data$experimental_condition, "_", combined_obj@meta.data$batch)

object = Seurat::NormalizeData(object)

wt_norm = object@assays$RNA@data[, rownames(object@meta.data[object@meta.data$experimental_condition == 'Wt', ])]
wt_mean_norm = apply(wt_norm, FUN = mean, MARGIN = 1)

mut_norm = object@assays$RNA@data[, rownames(object@meta.data[object@meta.data$experimental_condition == 'mut', ])]
mut_mean_norm = apply(mut_norm, FUN = mean, MARGIN = 1)

plot_df = data.frame(mut_norm = mut_mean_norm, 
                     wt_norm = wt_mean_norm)

down_genes = read.csv("../CrebA_ChIP_Analysis/output/find_bound_DE_genes/down_DE.csv", row.names = 1)
up_genes = read.csv("../CrebA_ChIP_Analysis/output/find_bound_DE_genes/up_DE.csv", row.names = 1)

down_genes = down_genes[down_genes$bound == 'True', ]
up_genes = up_genes[up_genes$bound == 'True', ]

down_genes = down_genes[down_genes$SC_DE == 'True', ]
up_genes = up_genes[up_genes$SC_DE == 'True', ]

plot_df$type = 'None'
plot_df[up_genes$genes, 'type'] = 'Repressed'
plot_df[down_genes$genes, 'type'] = 'Activated'

plot_df = plot_df[plot_df$mut_norm > 0 | plot_df$wt_norm > 0, ]

fit <- lm(mut_norm ~ wt_norm, data = plot_df)

# Extract slope and intercept
slope <- coef(fit)[2]
intercept <- coef(fit)[1]
r_squared <- summary(fit)$r.squared

# Create the regression equation as a label
equation <- paste0("y = ", round(slope, 2), "x + ", round(intercept, 2), 
                   ", R² = ", round(r_squared, 3))

p = ggplot(plot_df, aes(color = type)) +
  geom_point(data = subset(plot_df, type == "None"), aes(x = wt_norm, y = mut_norm), 
             color = "#ff7f0e", size = 3) +
  geom_point(data = subset(plot_df, type == "Activated"), aes(x = wt_norm, y = mut_norm), 
             color = "#1f77b4", size = 4) +
  geom_point(data = subset(plot_df, type == "Repressed"), aes(x = wt_norm, y = mut_norm), 
             color = "#2ca02c", size = 4) +
  geom_smooth(data = plot_df, aes(x = wt_norm, y = mut_norm), method = "lm", color = "black", linetype = "dotted", se = FALSE) +
  labs(
    title = "wt vs mut salivary gland norm experssion",
    x = "wt salivary gland norm expression",
    y = "mut salivary gland norm expression"
  ) +
  annotate("text", x = 0.5, 
           y = 4, 
           label = equation, 
           color = "black", size = 5, hjust = 0) +
  scale_color_manual(values = c("Activated" = "#1f77b4", "None" = "#ff7f0e", "Repressed" = "#2ca02c"), name = 'Legend') +
  cowplot::theme_cowplot()

lm_exp = lm(mut_norm ~ wt_norm, data = plot_df) #Create a linear regression with two variables
summary(lm_exp) #Review the results

ggsave(filename = file.path(TARGET_dir, 'wt_mut_norm_exp.png'), plot = p, width = 6, height = 5)

##### look at plasmatocyte #####

for(celltype in c('Salivary Gland', 'Plasmatocytes', 'Fat Body', 'Amnioserosa')) {
  sub_wt_obj = subset(wt_object, subset = manual_celltypes == celltype) 
  sub_wt_obj@meta.data$experimental_condition = 'Wt'
  
  sub_mut_obj = subset(mut_object, subset = manual_celltypes == celltype) 
  sub_mut_obj@meta.data$experimental_condition = 'mut'
  
  combined_obj = merge(sub_wt_obj, sub_mut_obj)
  object = Seurat::CreateSeuratObject(combined_obj@assays$RNA@counts, project = 'find_diff')
  object@meta.data$experimental_condition = combined_obj@meta.data$experimental_condition
  object@meta.data$batch = paste0(combined_obj@meta.data$experimental_condition, "_", combined_obj@meta.data$batch)
  
  object = Seurat::NormalizeData(object)
  
  wt_norm = object@assays$RNA@data[, rownames(object@meta.data[object@meta.data$experimental_condition == 'Wt', ])]
  wt_mean_norm = apply(wt_norm, FUN = mean, MARGIN = 1)
  
  mut_norm = object@assays$RNA@data[, rownames(object@meta.data[object@meta.data$experimental_condition == 'mut', ])]
  mut_mean_norm = apply(mut_norm, FUN = mean, MARGIN = 1)
  
  plot_df = data.frame(mut_norm = mut_mean_norm, 
                       wt_norm = wt_mean_norm)
  
  down_genes = read.csv("../CrebA_ChIP_Analysis/output/find_bound_DE_genes_majority_4cts/down_DE.csv", row.names = 1)
  up_genes = read.csv("../CrebA_ChIP_Analysis/output/find_bound_DE_genes_majority_4cts/up_DE.csv", row.names = 1)
  
  DE_genes = read.csv(file.path("results", ANALYSIS_VERSION, "/DE_genes_early_crebA_wt/", celltype,"/mut_DE_genes.csv"), row.names = 1)
  DE_genes = DE_genes[DE_genes$pval < 0.05 & abs(DE_genes$logFC) > 0.15, ]
  
  down_genes = down_genes[down_genes$genes %in% DE_genes$feature, ]
  up_genes = up_genes[up_genes$genes %in% DE_genes$feature, ]
  
  down_genes = down_genes[down_genes$bound == 'True', ]
  up_genes = up_genes[up_genes$bound == 'True', ]
  
  down_genes = down_genes[down_genes$SC_DE == 'True', ]
  up_genes = up_genes[up_genes$SC_DE == 'True', ]
  
  plot_df$type = 'None'
  plot_df[up_genes$genes, 'type'] = 'Repressed'
  plot_df[down_genes$genes, 'type'] = 'Activated'
  
  plot_df = plot_df[plot_df$mut_norm > 0 | plot_df$wt_norm > 0, ]
  
  fit <- lm(mut_norm ~ wt_norm, data = plot_df)
  
  # Extract slope and intercept
  slope <- coef(fit)[2]
  intercept <- coef(fit)[1]
  r_squared <- summary(fit)$r.squared
  
  # Create the regression equation as a label
  equation <- paste0("y = ", round(slope, 2), "x + ", round(intercept, 2), 
                     ", R² = ", round(r_squared, 3))
  
  p = ggplot(plot_df, aes(color = type)) +
    geom_point(data = subset(plot_df, type == "None"), aes(x = wt_norm, y = mut_norm), 
               color = "#ff7f0e", size = 3) +
    geom_point(data = subset(plot_df, type == "Activated"), aes(x = wt_norm, y = mut_norm), 
               color = "#1f77b4", size = 4) +
    geom_point(data = subset(plot_df, type == "Repressed"), aes(x = wt_norm, y = mut_norm), 
               color = "#2ca02c", size = 4) +
    geom_smooth(data = plot_df, aes(x = wt_norm, y = mut_norm), method = "lm", color = "black", linetype = "dotted", se = FALSE) +
    labs(
      title = paste0("wt vs mut ", celltype, " norm experssion"),
      x = paste0("wt ", celltype, " norm expression"),
      y = paste0("mut ", celltype, " norm expression")
    ) +
    annotate("text", x = 0.5, 
             y = 4, 
             label = equation, 
             color = "black", size = 5, hjust = 0) +
    scale_color_manual(values = c("Activated" = "#1f77b4", "None" = "#ff7f0e", "Repressed" = "#2ca02c"), name = 'Legend') +
    cowplot::theme_cowplot()
  
  lm_exp = lm(mut_norm ~ wt_norm, data = plot_df) #Create a linear regression with two variables
  summary(lm_exp) #Review the results
  
  ggsave(filename = file.path(TARGET_dir, paste0(celltype, '_wt_mut_norm_exp.png')), plot = p, width = 6, height = 5)
  
}

