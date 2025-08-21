TARGET_dir = file.path("results", ANALYSIS_VERSION, "Figures/plot_spcg_logFC")
dir.create(TARGET_dir)

##### load in the spcgs #####
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

# convert all the flybase name into the genes in seurat 
wt_object = readRDS(file.path("results", ANALYSIS_VERSION, "harmonized_wildtype_data", 'stage10-12_reharmonized_seurat.rds'))
flybase_map = rownames(wt_object)
spcg_tab$seurat_genes = NA
for(i in rownames(spcg_tab)) {
  spcg_tab[i, 'seurat_genes'] = flybase_map[spcg_tab[i, "Drosophila FBgn"]]
}
spcg_tab = spcg_tab[is.na(spcg_tab$seurat_genes) == FALSE, ]


new_row <- data.frame(
  "SPCG General Functional Categories" = "CrebA",
  "Drosophila Gene" = "CrebA",
  "Drosophila FBgn" = "CrebA",
  "seurat_genes" = 'CrebA', 
  check.names = FALSE
)

# Add the new row to the existing data frame
spcg_tab <- rbind(spcg_tab, new_row)
rownames(spcg_tab) = spcg_tab$seurat_genes

##### load and compile the early data #####
plot_df = data.frame()
ct_list = list.dirs(file.path('results', ANALYSIS_VERSION, 'DE_genes_early_crebA_wt'), full.names = FALSE, recursive = FALSE)
for(tmp_ct in ct_list) {
  mut_DE_genes = read.csv(file.path('results', ANALYSIS_VERSION, 'DE_genes_early_crebA_wt', tmp_ct, 'mut_DE_genes.csv'), row.names = 1)
  rownames(mut_DE_genes) = mut_DE_genes$feature
  mut_DE_genes = mut_DE_genes[mut_DE_genes$feature %in% spcg_tab$seurat_genes, ]
  mut_DE_genes$celltype = tmp_ct
  mut_DE_genes$spcg_cat = spcg_tab[mut_DE_genes$feature, "SPCG General Functional Categories"]
  plot_df = rbind(plot_df, mut_DE_genes)
}

exp_order = read.csv(file.path(file.path("results", ANALYSIS_VERSION, "Figures/plot_spcg_wt", 'wt_early', 'spcg_scale_exp.csv')), row.names = 1)
sub_exp_order = exp_order[exp_order$features.plot == 'CrebA', ]
sub_exp_order = sub_exp_order[sub_exp_order$id %in% plot_df$celltype, ]
sub_exp_order = sub_exp_order[order(sub_exp_order$avg.exp.scaled, decreasing = FALSE), ]

plot_df$celltype = factor(plot_df$celltype, levels = sub_exp_order$id)
plot_df$spcg_cat = factor(plot_df$spcg_cat, levels = c("CrebA", unique(spcg_tab$`SPCG General Functional Categories`)[unique(spcg_tab$`SPCG General Functional Categories`) != 'CrebA']))
plot_df = plot_df[plot_df$spcg_cat != 'Prolyl hydroxylation', ]

# order the genes 
sg_exp_order = exp_order[exp_order$id == 'Salivary Gland', ]
plot_df$feature = factor(plot_df$feature, levels = sg_exp_order$features.plot)
write.csv(plot_df, file = file.path(TARGET_dir, 'stage10-12_logFC.csv'))

p <- ggplot(data = plot_df, mapping = aes_string(y = 'celltype', x = 'feature', fill = 'logFC')) +
  geom_tile() +
  guides(fill = guide_colorbar(title = 'logFC', barwidth=30)) +
  scale_fill_gradient2(low = scales::muted("red"), mid = "white", high = scales::muted("blue"), 
                        midpoint = 0, limits = c(-1.6, 0.5)) +
  labs(
    x = 'Secretory Pathway Component Genes',
    y = 'Cell Types'
  ) + 
  facet_grid(
    cols = vars(spcg_cat),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=45, vjust = 1, hjust=1), 
        axis.text = element_text(size=18), text = element_text(size = 18), 
        legend.position = 'bottom') +
  ggtitle("Stage 10-12 Embryos")
ggsave(filename = file.path(TARGET_dir, "stage10-12_logFC.png"), plot = p, width = 32, height = 7)

# yellow and purple 
p <- ggplot(data = plot_df, mapping = aes_string(y = 'celltype', x = 'feature', fill = 'logFC')) +
  geom_tile() +
  guides(fill = guide_colorbar(title = 'logFC', barwidth=30)) +
  scale_fill_gradient2(low = "#d4ad00", mid = "white", high = '#6f02b8', 
                       midpoint = 0, limits = c(-1.6, 0.5)) +
  labs(
    x = 'Secretory Pathway Component Genes',
    y = 'Cell Types'
  ) + 
  facet_grid(
    cols = vars(spcg_cat),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=45, vjust = 1, hjust=1), 
        axis.text = element_text(size=18), text = element_text(size = 18), 
        legend.position = 'bottom') +
  ggtitle("Stage 10-12 Embryos")
ggsave(filename = file.path(TARGET_dir, "stage10-12_logFC_yellow_purple.png"), plot = p, width = 32, height = 7)

# purple and yellow
p <- ggplot(data = plot_df, mapping = aes_string(y = 'celltype', x = 'feature', fill = 'logFC')) +
  geom_tile() +
  guides(fill = guide_colorbar(title = 'logFC', barwidth=30)) +
  scale_fill_gradient2(low = "#6f02b8", mid = "white", high = '#d4ad00', 
                       midpoint = 0, limits = c(-1.6, 0.5)) +
  labs(
    x = 'Secretory Pathway Component Genes',
    y = 'Cell Types'
  ) + 
  facet_grid(
    cols = vars(spcg_cat),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=45, vjust = 1, hjust=1), 
        axis.text = element_text(size=18), text = element_text(size = 18), 
        legend.position = 'bottom') +
  ggtitle("Stage 10-12 Embryos")
ggsave(filename = file.path(TARGET_dir, "stage10-12_logFC_purple_yellow.png"), plot = p, width = 32, height = 7)

##### plot out the SPCG as box plot ######
color_palette = readRDS(file.path("results", ANALYSIS_VERSION, "ct_color_palettes", 'ct_color_palette.rds'))
spcg_plot_df = plot_df[plot_df$spcg_cat != 'CrebA', ]
spcg_plot_df$celltype = factor(spcg_plot_df$celltype, levels = rev(sub_exp_order$id))

p = ggplot(spcg_plot_df, aes(x = celltype, y = logFC, fill = celltype)) + 
  geom_boxplot() + 
  scale_fill_manual(values = color_palette, name = 'Cell Types') +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab('logFC of SPCGs') + 
  xlab('Cell Types')
ggsave(filename = file.path(TARGET_dir, 'stage10-12_logFC_boxplot.png'), plot = p, width = 12, height = 7)

##### load and compile the late data #####
plot_df = data.frame()
ct_list = list.dirs(file.path('results', ANALYSIS_VERSION, 'DE_genes_crebA_wt'), full.names = FALSE, recursive = FALSE)
for(tmp_ct in ct_list) {
  mut_DE_genes = read.csv(file.path('results', ANALYSIS_VERSION, 'DE_genes_crebA_wt', tmp_ct, 'mut_DE_genes.csv'), row.names = 1)
  rownames(mut_DE_genes) = mut_DE_genes$feature
  mut_DE_genes = mut_DE_genes[mut_DE_genes$feature %in% spcg_tab$seurat_genes, ]
  mut_DE_genes$celltype = tmp_ct
  mut_DE_genes$spcg_cat = spcg_tab[mut_DE_genes$feature, "SPCG General Functional Categories"]
  plot_df = rbind(plot_df, mut_DE_genes)
}

sub_exp_order = plot_df[plot_df$feature == 'CrebA', ]
sub_exp_order = sub_exp_order[order(sub_exp_order$logFC, decreasing = TRUE), ]

plot_df$celltype = factor(plot_df$celltype, levels = sub_exp_order$celltype)
plot_df$spcg_cat = factor(plot_df$spcg_cat, levels = c("CrebA", unique(spcg_tab$`SPCG General Functional Categories`)[unique(spcg_tab$`SPCG General Functional Categories`) != 'CrebA']))
plot_df = plot_df[plot_df$spcg_cat != 'Prolyl hydroxylation', ]

write.csv(plot_df, file = file.path(TARGET_dir, 'stage13-16_logFC.csv'))
p <- ggplot(data = plot_df, mapping = aes_string(y = 'celltype', x = 'feature', fill = 'logFC')) +
  geom_tile() +
  guides(fill = guide_colorbar(title = 'logFC', barwidth=30)) +
  scale_fill_gradient2(low = scales::muted("red"), mid = "white", high = scales::muted("blue"), 
                       midpoint = 0, limits = c(min(plot_df$logFC), max(plot_df$logFC))) +
  labs(
    x = 'Secretory Pathway Component Genes',
    y = 'Cell Types'
  ) + 
  facet_grid(
    cols = vars(spcg_cat),
    scales = "free_x",
    space = "free_x",
    switch = "y"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank()
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=45, vjust = 1, hjust=1), 
        axis.text = element_text(size=18), text = element_text(size = 18), 
        legend.position = 'bottom') +
  ggtitle("Stage 13-16 Embryos")
ggsave(filename = file.path(TARGET_dir, "stage13-16_logFC.png"), plot = p, width = 31, height = 7)

