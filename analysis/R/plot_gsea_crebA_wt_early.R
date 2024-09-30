TARGET_dir = file.path("results", ANALYSIS_VERSION, "Figures/plot_gsea_early_crebA_wt")
dir.create(TARGET_dir)

##### write a function to automatically get the top 5 categories for each cell type #####
get_top_cat <- function(ct, top_num = 5) {
  wt_gsea = read.csv(file.path("results", ANALYSIS_VERSION, "DE_genes_early_crebA_wt", ct, 'wt_gsea_results.csv'), row.names = 1)
  mut_gsea = read.csv(file.path("results", ANALYSIS_VERSION, "DE_genes_early_crebA_wt", ct, 'mut_gsea_results.csv'), row.names = 1)
  
  wt_gsea = wt_gsea[wt_gsea$NES > 0, ]
  mut_gsea = mut_gsea[mut_gsea$NES > 0, ]
  
  wt_gsea = wt_gsea[order(wt_gsea$NES, decreasing = TRUE), ]
  mut_gsea = mut_gsea[order(mut_gsea$NES, decreasing = TRUE), ]
  
  wt_gsea = wt_gsea[1:top_num, ]
  mut_gsea = mut_gsea[1:top_num, ]
  
  wt_gsea$type = 'wildtype'
  mut_gsea$type = 'mutant'
  
  wt_gsea$logpval = -log(wt_gsea$padj)
  mut_gsea$logpval = -log(wt_gsea$padj)
  return(rbind(wt_gsea, mut_gsea))
}

ct_list = list.dirs(file.path("results", ANALYSIS_VERSION, 'DE_genes_early_crebA_wt'), full.names = FALSE, recursive = FALSE)
for(ct in ct_list) {
  plot_df = get_top_cat(ct)
  p <- ggplot(data = plot_df, aes(y = reorder(pathway, logpval), x = logpval, fill = type)) +
    geom_bar(stat="identity") +
    labs(
      x = '-log10 adjusted p-value',
      y = ''
    ) + 
    theme_classic()  + 
    facet_grid(
      rows = vars(type),
      scales = "free_y",
      space = "free_y",
      switch = "x"
    ) + 
    theme(
      panel.spacing = unit(x = 1, units = "lines"),
      strip.background = element_blank(), 
      text = element_text(size = 25), 
      legend.position="none"
    ) + 
    theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=0, vjust = 1, hjust=1)) +
    ggtitle(ct)  
  ggsave(filename = file.path(TARGET_dir, paste0(ct, "_gsea.png")), width = 15, height = 8)
}

##### write the function that curate all the data ######
compiled_df = data.frame()
for(ct in ct_list) {
  tmp_df = get_top_cat(ct, top_num = 4)
  tmp_df$celltype = ct
  compiled_df = rbind(compiled_df, tmp_df)
}

wt_df = compiled_df[compiled_df$type == 'wildtype', ]
mut_df = compiled_df[compiled_df$type == 'mutant', ]

##### make the wildtype plot #####
color_palette = readRDS(file.path('results', ANALYSIS_VERSION, 'ct_color_palettes/ct_color_palette.rds'))
sub_wt_df = wt_df[wt_df$celltype %in% c('Salivary Gland', 'Plasmatocytes', 'Amnioserosa', 'Fat Body'), ]
sub_wt_df$celltype = factor(x = sub_wt_df$celltype, levels = c('Salivary Gland', 'Amnioserosa', 'Plasmatocytes', 'Fat Body'))
sub_wt_df$pathway = stringr::str_split_fixed(sub_wt_df$pathway, pattern = " \\(", n = 2)[, 1]

p <- ggplot(data = sub_wt_df, aes(y = reorder(pathway, logpval), x = logpval, fill = celltype)) +
  geom_bar(stat="identity") +
  labs(
    x = '-log10 adjusted p-value',
    y = ''
  ) + 
  scale_fill_manual(values =color_palette) + 
  theme_classic()  + 
  facet_grid(
    rows = vars(celltype),
    scales = "free_y",
    space = "free_y",
    switch = "x"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank(), 
    text = element_text(size = 25), 
    legend.position="none", 
    plot.title.position = "plot"
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=0, vjust = 1, hjust=1)) +
  ggtitle('Geneset enrichment in stage 10-12 wildtype embryos')  
ggsave(filename = file.path(TARGET_dir, 'wildtype_selected_ct.png'), width = 11, height = 10)

sub_mut_df = mut_df[mut_df$celltype %in% c('Salivary Gland', 'Plasmatocytes', 'Amnioserosa', 'Fat Body'), ]
sub_mut_df$celltype = factor(x = sub_mut_df$celltype, levels = c('Salivary Gland', 'Amnioserosa', 'Plasmatocytes', 'Fat Body'))

sub_mut_df$pathway = stringr::str_split_fixed(sub_mut_df$pathway, pattern = " \\(", n = 2)[, 1]

p <- ggplot(data = sub_mut_df, aes(y = reorder(pathway, logpval), x = logpval, fill = celltype)) +
  geom_bar(stat="identity") +
  labs(
    x = '-log10 adjusted p-value',
    y = ''
  ) + 
  scale_fill_manual(values =color_palette) + 
  theme_classic()  + 
  facet_grid(
    rows = vars(celltype),
    scales = "free_y",
    space = "free_y",
    switch = "x"
  ) + 
  theme(
    panel.spacing = unit(x = 1, units = "lines"),
    strip.background = element_blank(), 
    text = element_text(size = 25), 
    legend.position="none", 
    plot.title.position = "plot"
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=0, vjust = 1, hjust=1)) +
  ggtitle('Geneset enrichment in stage 10-12 mutant embryos')  
ggsave(filename = file.path(TARGET_dir, 'mutant_selected_ct.png'), width = 11, height = 10)

