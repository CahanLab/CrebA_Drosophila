TARGET_dir = file.path("results", ANALYSIS_VERSION, "Figures/plot_spcg_wt")
dir.create(TARGET_dir)

spcg_tab = readxl::read_excel("accessory_data/SPCG_files/SPCG List.xlsx")
spcg_tab = as.data.frame(spcg_tab)
for(i in seq(1, nrow(spcg_tab))) {
  if(is.na(spcg_tab[i, 'SPCG General Functional Categories']) == TRUE) {
    spcg_tab[i, 'SPCG General Functional Categories'] = spcg_tab[i - 1, 'SPCG General Functional Categories'] 
  }
}

spcg_tab = spcg_tab[spcg_tab$`Drosophila Gene` %in% c('CG9356',
                                                      'Cog1',
                                                      'Cog2',
                                                      'Cog3',
                                                      'Cog6',
                                                      'Cog7',
                                                      'fws',
                                                      'Osbp',
                                                      'p24-2',
                                                      'PH4alphaSG1',
                                                      'rt',
                                                      'Sec71',
                                                      'Sil1', 
                                                      'spas',
                                                      'TBC1D23', 
                                                      'tw') == FALSE, ]
##### load in the color palette ######
color_palette = readRDS(file.path("results", ANALYSIS_VERSION, "ct_color_palettes", 'ct_color_palette.rds'))

##### make the early wildtype plots for the SPCG genes ######
wt_object = readRDS(file.path("results", ANALYSIS_VERSION, "harmonized_wildtype_data", 'stage10-12_reharmonized_seurat.rds'))
subfolder_path = 'wt_early'
if(dir.exists(file.path(TARGET_dir, subfolder_path)) == FALSE) {
  dir.create(file.path(TARGET_dir, subfolder_path))
}

flybase_map = rownames(wt_object)
spcg_tab$seurat_genes = NA
for(i in rownames(spcg_tab)) {
  spcg_tab[i, 'seurat_genes'] = flybase_map[spcg_tab[i, "Drosophila FBgn"]]
}
spcg_tab = spcg_tab[is.na(spcg_tab$seurat_genes) == FALSE, ]

# this is to reorder violin plots based on CrebA expression
mean_expression = AverageExpression(wt_object, features = "CrebA", group.by = 'new_celltypes')$RNA[1, ]
ordered_clusters = order(mean_expression, decreasing = TRUE)
ordered_cluster_names = names(mean_expression)[ordered_clusters]

# Create a new factor level with the ordered clusters
wt_object$new_celltypes = factor(wt_object$new_celltypes, levels = ordered_cluster_names)

# this is to plot the violin plots 
for(gene_id in c(spcg_tab$seurat_genes, 'CrebA')) {
  p = VlnPlot(wt_object, features = gene_id, group.by = 'new_celltypes', pt.size = 0) + 
      scale_fill_manual(values = color_palette) +
      ggtitle(paste0("Stage 10-12 wildtype embryos: ", gene_id, " expression")) + 
      theme(legend.position = "none", axis.text = element_text(size=17), text = element_text(size = 17), plot.margin = margin(0, 0, 0, 1.6, "cm")) +
      geom_boxplot(width=0.2, color="black", fill = 'white') + 
      xlab("Cell Types") + 
      ylab("Expression Level")
    ggsave(filename = file.path(TARGET_dir, subfolder_path, 'vlnPlots', paste0(gene_id, "_early_vlnplot.png")), plot = p, width = 16, height = 6)
}

# get the genes categorized 
spcg_dot_df = modified_dotPlot_df(wt_object, features = unique(c(spcg_tab$seurat_genes, 'CrebA')), group.by = 'new_celltypes')
spcg_dot_df$category = NA

for(temp_gene in unique(spcg_dot_df$features.plot)) { 
  if(temp_gene == 'CrebA') {
    cat = 'CrebA'
  } else { 
    sub_spcg_tab = spcg_tab[spcg_tab$seurat_genes == temp_gene, ]
    cat = paste(sub_spcg_tab[, "SPCG General Functional Categories"], collapse = ', ')
  }
  spcg_dot_df[spcg_dot_df$features.plot == temp_gene, 'category'] = cat
  
}

sub_crebA_exp = spcg_dot_df[spcg_dot_df$features.plot == 'CrebA', ]
sub_crebA_exp = sub_crebA_exp[order(sub_crebA_exp$avg.exp.scaled), ]

spcg_dot_df$id = factor(spcg_dot_df$id, levels = sub_crebA_exp$id)
spcg_dot_df$category = factor(spcg_dot_df$category, levels = c('CrebA', unique(spcg_dot_df$category)[unique(spcg_dot_df$category) != 'CrebA']))
spcg_dot_df = spcg_dot_df[spcg_dot_df$category != 'Prolyl hydroxylation', ]

write.csv(spcg_dot_df, file = file.path(TARGET_dir, subfolder_path, 'spcg_scale_exp.csv'))

# plot out the dot plot 
p <- ggplot(data = spcg_dot_df, mapping = aes_string(y = 'id', x = 'features.plot')) +
  geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp.scaled')) +
  guides(size = guide_legend(title = 'percent expressed')) +
  guides(color = guide_colorbar(title = 'scaled average expression')) +
  scale_colour_viridis_c() + 
  labs(
    x = 'Secretory Pathway Component Genes',
    y = 'Cell Types'
  ) + 
  theme_classic()  + 
  facet_grid(
    cols = vars(category),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.x = element_blank(), 
        axis.text.x=element_text(angle=45, vjust = 1, hjust=1), 
        axis.text = element_text(size=17), text = element_text(size = 17), 
        legend.position = 'bottom') +
  ggtitle("Stage 10-12 Embryos")
ggsave(filename = file.path(TARGET_dir, subfolder_path, "spcg_dot_scaled.png"), plot = p, width = 32, height = 10)

# plot out the unscaled version 
p <- ggplot(data = spcg_dot_df, mapping = aes_string(y = 'id', x = 'features.plot')) +
  geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) +
  guides(size = guide_legend(title = 'percent expressed')) +
  guides(color = guide_colorbar(title = 'average expression')) +
  scale_colour_viridis_c() + 
  labs(
    x = 'Secretory Pathway Component Genes',
    y = 'Cell Types'
  ) + 
  theme_classic()  + 
  facet_grid(
    cols = vars(category),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=45, vjust = 1, hjust=1), 
        axis.text = element_text(size=15), text = element_text(size = 15), 
        legend.position = 'bottom') +
  ggtitle("Stage 10-12 Embryos")
ggsave(filename = file.path(TARGET_dir, subfolder_path, "spcg_dot_norm_exp.png"), plot = p, width = 31, height = 7)

##### this is for late wild type samples ######
wt_object = readRDS(file.path("results", ANALYSIS_VERSION, "harmonized_wildtype_data", 'stage13-16_reharmonized_seurat.rds'))
subfolder_path = 'wt_late'
if(dir.exists(file.path(TARGET_dir, subfolder_path)) == FALSE) {
  dir.create(file.path(TARGET_dir, subfolder_path))
}

flybase_map = rownames(wt_object)
spcg_tab$seurat_genes = NA
for(i in rownames(spcg_tab)) {
  spcg_tab[i, 'seurat_genes'] = flybase_map[spcg_tab[i, "Drosophila FBgn"]]
}
spcg_tab = spcg_tab[is.na(spcg_tab$seurat_genes) == FALSE, ]

# this is to reorder violin plots based on CrebA expression
mean_expression = AverageExpression(wt_object, features = "CrebA", group.by = 'new_celltypes')$RNA[1, ]
ordered_clusters = order(mean_expression, decreasing = TRUE)
ordered_cluster_names = names(mean_expression)[ordered_clusters]

# Create a new factor level with the ordered clusters
wt_object$new_celltypes = factor(wt_object$new_celltypes, levels = ordered_cluster_names)

# this is to plot the violin plots 
for(gene_id in c(spcg_tab$seurat_genes, 'CrebA')) {
  p = VlnPlot(wt_object, features = gene_id, group.by = 'new_celltypes', pt.size = 0) + 
    scale_fill_manual(values = color_palette) +
    ggtitle(paste0("Stage 13-16 wildtype embryos: ", gene_id, " expression")) + 
    theme(legend.position = "none", axis.text = element_text(size=17), text = element_text(size = 17), plot.margin = margin(0, 0, 0, 1.6, "cm")) +
    geom_boxplot(width=0.2, color="black", fill = 'white') + 
    xlab("Cell Types") + 
    ylab("Expression Level")
  ggsave(filename = file.path(TARGET_dir, subfolder_path, 'vlnPlots', paste0(gene_id, "_late_vlnplot.png")), plot = p, width = 18, height = 6)
}

# get the genes categorized 
spcg_dot_df = modified_dotPlot_df(wt_object, features = unique(c(spcg_tab$seurat_genes, 'CrebA')), group.by = 'manual_celltypes')
spcg_dot_df$category = NA

for(temp_gene in unique(spcg_dot_df$features.plot)) { 
  if(temp_gene == 'CrebA') {
    cat = 'CrebA'
  } else { 
    sub_spcg_tab = spcg_tab[spcg_tab$seurat_genes == temp_gene, ]
    cat = paste(sub_spcg_tab[, "SPCG General Functional Categories"], collapse = ', ')
  }
  spcg_dot_df[spcg_dot_df$features.plot == temp_gene, 'category'] = cat
  
}
sub_crebA_exp = spcg_dot_df[spcg_dot_df$features.plot == 'CrebA', ]
sub_crebA_exp = sub_crebA_exp[order(sub_crebA_exp$avg.exp.scaled), ]

spcg_dot_df$id = factor(spcg_dot_df$id, levels = sub_crebA_exp$id)
spcg_dot_df$category = factor(spcg_dot_df$category, levels = unique(spcg_dot_df$category))
spcg_dot_df = spcg_dot_df[spcg_dot_df$category != 'Prolyl hydroxylation', ]

write.csv(spcg_dot_df, file = file.path(TARGET_dir, subfolder_path, 'spcg_scale_exp.csv'))

# plot out the dot plot 
p <- ggplot(data = spcg_dot_df, mapping = aes_string(y = 'id', x = 'features.plot')) +
  geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp.scaled')) +
  guides(size = guide_legend(title = 'percent expressed')) +
  guides(color = guide_colorbar(title = 'scaled average expression')) +
  scale_colour_viridis_c() + 
  labs(
    x = 'Secretory Pathway Component Genes',
    y = 'Cell Types'
  ) + 
  theme_classic()  + 
  facet_grid(
    cols = vars(category),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=45, vjust = 1, hjust=1), 
        axis.text = element_text(size=15), text = element_text(size = 15), 
        legend.position = 'bottom') +
  ggtitle("Stage 13-16 Embryos")
ggsave(filename = file.path(TARGET_dir, subfolder_path, "spcg_dot_scaled.png"), plot = p, width = 32, height = 8)

# plot out the dot plot 
p <- ggplot(data = spcg_dot_df, mapping = aes_string(y = 'id', x = 'features.plot')) +
  geom_point(mapping = aes_string(size = 'pct.exp', color = 'avg.exp')) +
  guides(size = guide_legend(title = 'percent expressed')) +
  guides(color = guide_colorbar(title = 'average expression')) +
  scale_colour_viridis_c() + 
  labs(
    x = 'Secretory Pathway Component Genes',
    y = 'Cell Types'
  ) + 
  theme_classic()  + 
  facet_grid(
    cols = vars(category),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=45, vjust = 1, hjust=1), 
        axis.text = element_text(size=15), text = element_text(size = 15), 
        legend.position = 'bottom') +
  ggtitle("Stage 13-16 Embryos")
ggsave(filename = file.path(TARGET_dir, subfolder_path, "spcg_dot_norm_exp.png"), plot = p, width = 32, height = 8)
