TARGET_dir = file.path("results", ANALYSIS_VERSION, "Figures/plot_secretome_wt")
dir.create(TARGET_dir)

secretome_tab = readxl::read_excel("accessory_data/curated_secretome/WT Salivary gland_Secretome.xlsx", sheet = 'Secretome list-functions')
secretome_tab = as.data.frame(secretome_tab)

##### load in the color palette ######
color_palette = readRDS(file.path("results", ANALYSIS_VERSION, "ct_color_palettes", 'ct_color_palette.rds'))

##### make the early wildtype plots for the SPCG genes ######
wt_object = readRDS(file.path("results", ANALYSIS_VERSION, "harmonized_wildtype_data", 'stage10-12_reharmonized_seurat.rds'))
subfolder_path = 'wt_early'
if(dir.exists(file.path(TARGET_dir, subfolder_path)) == FALSE) {
  dir.create(file.path(TARGET_dir, subfolder_path))
}

flybase_map = rownames(wt_object)
secretome_tab$seurat_genes = NA
for(i in rownames(secretome_tab)) {
  secretome_tab[i, 'seurat_genes'] = flybase_map[secretome_tab[i, "flybase_id"]]
}
secretome_tab = secretome_tab[is.na(secretome_tab$seurat_genes) == FALSE, ]

presto_results = presto::wilcoxauc(wt_object, group_by = 'new_celltypes')
presto_results = presto_results[presto_results$feature %in% secretome_tab$seurat_genes, ]
presto_results = presto_results[presto_results$pct_in > 10, ]
secretome_tab = secretome_tab[secretome_tab$seurat_genes %in% presto_results$feature, ]

# this is to plot the violin plots 
for(gene_id in secretome_tab$seurat_genes) {
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
secretome_dot_df = modified_dotPlot_df(wt_object, features = unique(secretome_tab$seurat_genes), group.by = 'new_celltypes')
secretome_dot_df$category = NA

for(temp_gene in unique(secretome_dot_df$features.plot)) { 

  sub_secretome_tab = secretome_tab[secretome_tab$seurat_genes == temp_gene, ]
  cat = paste(sub_secretome_tab[, "Biological or molecular function"], collapse = ', ')

  secretome_dot_df[secretome_dot_df$features.plot == temp_gene, 'category'] = cat
  
}

write.csv(secretome_dot_df, file = file.path(TARGET_dir, subfolder_path, 'spcg_scale_exp.csv'))

# plot out the dot plot 
p <- ggplot(data = secretome_dot_df, mapping = aes_string(y = 'id', x = 'features.plot')) +
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
ggsave(filename = file.path(TARGET_dir, subfolder_path, "spcg_dot_scaled.png"), plot = p, width = 49, height = 10)

# plot out the unscaled version 
p <- ggplot(data = secretome_dot_df, mapping = aes_string(y = 'id', x = 'features.plot')) +
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
ggsave(filename = file.path(TARGET_dir, subfolder_path, "spcg_dot_norm_exp.png"), plot = p, width = 49, height = 10)

##### this is for late wild type samples ######
wt_object = readRDS(file.path("results", ANALYSIS_VERSION, "harmonized_wildtype_data", 'stage13-16_reharmonized_seurat.rds'))
subfolder_path = 'wt_late'
if(dir.exists(file.path(TARGET_dir, subfolder_path)) == FALSE) {
  dir.create(file.path(TARGET_dir, subfolder_path))
}

flybase_map = rownames(wt_object)
secretome_tab$seurat_genes = NA
for(i in rownames(secretome_tab)) {
  secretome_tab[i, 'seurat_genes'] = flybase_map[secretome_tab[i, "flybase_id"]]
}
secretome_tab = secretome_tab[is.na(secretome_tab$seurat_genes) == FALSE, ]

presto_results = presto::wilcoxauc(wt_object, group_by = 'new_celltypes')
presto_results = presto_results[presto_results$feature %in% secretome_tab$seurat_genes, ]
presto_results = presto_results[presto_results$pct_in > 10, ]
secretome_tab = secretome_tab[secretome_tab$seurat_genes %in% presto_results$feature, ]

# this is to plot the violin plots 
for(gene_id in secretome_tab$seurat_genes) {
  p = VlnPlot(wt_object, features = gene_id, group.by = 'new_celltypes', pt.size = 0) + 
    scale_fill_manual(values = color_palette) +
    ggtitle(paste0("Stage 13-16 wildtype embryos: ", gene_id, " expression")) + 
    theme(legend.position = "none", axis.text = element_text(size=17), text = element_text(size = 17), plot.margin = margin(0, 0, 0, 1.6, "cm")) +
    geom_boxplot(width=0.2, color="black", fill = 'white') + 
    xlab("Cell Types") + 
    ylab("Expression Level")
  ggsave(filename = file.path(TARGET_dir, subfolder_path, 'vlnPlots', paste0(gene_id, "_late_vlnplot.png")), plot = p, width = 16, height = 6)
}

# get the genes categorized 
secretome_dot_df = modified_dotPlot_df(wt_object, features = unique(secretome_tab$seurat_genes), group.by = 'new_celltypes')
secretome_dot_df$category = NA

for(temp_gene in unique(secretome_dot_df$features.plot)) { 
  
  sub_secretome_tab = secretome_tab[secretome_tab$seurat_genes == temp_gene, ]
  cat = paste(sub_secretome_tab[, "Biological or molecular function"], collapse = ', ')
  
  secretome_dot_df[secretome_dot_df$features.plot == temp_gene, 'category'] = cat
  
}

write.csv(secretome_dot_df, file = file.path(TARGET_dir, subfolder_path, 'spcg_scale_exp.csv'))

# plot out the dot plot 
p <- ggplot(data = secretome_dot_df, mapping = aes_string(y = 'id', x = 'features.plot')) +
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
  ggtitle("Stage 13-16 Embryos")
ggsave(filename = file.path(TARGET_dir, subfolder_path, "spcg_dot_scaled.png"), plot = p, width = 49, height = 10)

# plot out the unscaled version 
p <- ggplot(data = secretome_dot_df, mapping = aes_string(y = 'id', x = 'features.plot')) +
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
ggsave(filename = file.path(TARGET_dir, subfolder_path, "spcg_dot_norm_exp.png"), plot = p, width = 49, height = 10)
