TARGET_dir = file.path("results", ANALYSIS_VERSION, "categorize_DE_genes")
dir.create(TARGET_dir)

###### identify the cell types with sufficient CrebA #####
wt_object_1 = readRDS("accessory_data/wild_type_seurats/stage_10-12/manual_celltype_object1.rds")
wt_object_2 = readRDS("accessory_data/wild_type_seurats/stage_13-16/manual_celltype_object4.rds")

# select the early wildtype 
meta_object_1 = wt_object_1@meta.data
meta_object_1$CrebA = wt_object_1@assays$RNA@data['CrebA', ]
meta_median_1 = meta_object_1 %>%
                group_by(manual_celltypes) %>% 
                summarise(median_values = median(CrebA, na.rm = TRUE))
meta_median_1 = meta_median_1[meta_median_1$median_values > 1, ]

# select the late wildtype 
meta_object_2 = wt_object_2@meta.data
meta_object_2$CrebA = wt_object_2@assays$RNA@data['CrebA', ]
meta_median_2 = meta_object_2 %>% 
  group_by(manual_celltypes) %>% 
  summarise(median_values = median(CrebA, na.rm = TRUE))
meta_median_2 = meta_median_2[meta_median_2$median_values > 1, ]

###### load in the microarray data #####
MA_data = readxl::read_excel("accessory_data/CrebA_microarry_data/CrebA_all microarray data.xlsx")
MA_data = MA_data[MA_data$`Gene Symbol` != '---', ]
MA_data$`Gene Symbol` = stringr::str_split_fixed(MA_data$`Gene Symbol`, " ///", n = 2) %>% `[`(, 1)
MA_data = MA_data[!duplicated(MA_data$`Gene Symbol`), ]
MA_data$logFC = log2(abs(MA_data$`CrebA vs. WT linear FC`)) * sign(MA_data$`CrebA vs. WT linear FC`)
MA_data = data.frame(MA_data)
rownames(MA_data) = MA_data$Gene.Symbol

###### load in the CrebA ChIP-data that match with the closest gene ######
intersect_ChIP_match = read.csv("accessory_data/CrebA_ChIP_Analysis/output/match_nearest_intersect_250_peaks/match_nearest_genes.csv", row.names = 1)
wildtype_ChIP_match = read.csv("accessory_data/CrebA_ChIP_Analysis/output/match_nearest_wildtype_250_peaks/match_nearest_genes.csv", row.names = 1)

###### get the scRNA logFC for early ######
save_path = file.path(TARGET_dir, 'early_data')
dir.create(save_path)
dir.create(file.path(save_path, 'plots_ct'))
dir.create(file.path(save_path, 'selected_genes_ct'))
full_gene_df = data.frame()
for(tmp_cat in unique(meta_median_1$manual_celltypes)) {
  print(tmp_cat)
  SC_data = read.csv(file.path("results/v19/DE_genes_early_crebA_wt/", tmp_cat, 'mut_DE_genes.csv'), row.names = 1)
  rownames(SC_data) = SC_data$feature
  sub_MA_data = MA_data[intersect(SC_data$feature, MA_data$Gene.Symbol), ]
  sub_SC_data = SC_data[intersect(SC_data$feature, MA_data$Gene.Symbol), ]
  tmp_df = data.frame('ct' = tmp_cat, 
                      'gene' = intersect(SC_data$feature, MA_data$Gene.Symbol),
                      'sc_logFC' = sub_SC_data$logFC, 
                      'ma_logFC' = sub_MA_data$logFC)
  full_gene_df = rbind(full_gene_df, tmp_df)
}

# look at the down up categories across different cell types 
interesting_genes = c()
for(tmp_ct in unique(full_gene_df$ct)) {
  subset_gene_df = full_gene_df[full_gene_df$ct == tmp_ct, ]
  subset_gene_df$categories = 'other'
  subset_gene_df[subset_gene_df$sc_logFC < -0.2 & subset_gene_df$ma_logFC < -0.3, 'categories'] = 'down'
  subset_gene_df[subset_gene_df$sc_logFC > 0.2 & subset_gene_df$ma_logFC > 0.3, 'categories'] = 'up'
  subset_gene_df[subset_gene_df$sc_logFC > -0.05 & subset_gene_df$sc_logFC < 0.05 &
                   subset_gene_df$ma_logFC > -0.05 & subset_gene_df$ma_logFC < 0.05, 'categories'] = 'static'
  p = ggplot(data = subset_gene_df, aes(x = ma_logFC, y = sc_logFC, color = categories)) +
    geom_point() +
    theme_cowplot() +
    labs(title = "Early single-cell vs MA",
         x = "MA DE Genes LogFC",
         y = "SC DE Genes logFC") + 
    ggtitle(tmp_ct)
  ggsave(file.path(save_path, 'plots_ct', paste0(tmp_ct, '.png')), plot = p, width = 7, height = 5)
  subset_gene_df = subset_gene_df[subset_gene_df$categories != 'other', ]
  interesting_genes = c(interesting_genes, subset_gene_df$gene)
  write.csv(subset_gene_df, file = file.path(save_path, 'selected_genes_ct', paste0(tmp_ct, '.csv')))
}

# get the status table across the cell types 
interesting_genes = unique(interesting_genes)
cat_df = data.frame(genes = interesting_genes)
for(tmp_ct in unique(full_gene_df$ct)) {
  subset_gene_df = full_gene_df[full_gene_df$ct == tmp_ct, ]
  subset_gene_df$categories = 'other'
  subset_gene_df[subset_gene_df$sc_logFC < -0.2 & subset_gene_df$ma_logFC < -0.3, 'categories'] = 'down'
  subset_gene_df[subset_gene_df$sc_logFC > 0.2 & subset_gene_df$ma_logFC > 0.3, 'categories'] = 'up'
  subset_gene_df[subset_gene_df$sc_logFC > -0.05 & subset_gene_df$sc_logFC < 0.05 &
                   subset_gene_df$ma_logFC > -0.05 & subset_gene_df$ma_logFC < 0.05, 'categories'] = 'static'
  rownames(subset_gene_df) = subset_gene_df$gene
  cat_df[tmp_ct] = subset_gene_df[interesting_genes, 'categories']
}
cat_df[is.na(cat_df)] = 'other'

# check if the categorized genes have CrebA binding either salivary gland specific or overall
cat_df$fkh_sage_CrebA_bind = 0
cat_df[cat_df$genes %in% c(intersect_ChIP_match$nearest_gene, intersect_ChIP_match$nearest_gene2), 'fkh_sage_CrebA_bind'] = 1
cat_df$wt_CrebA_bind = 0
cat_df[cat_df$genes %in% c(wildtype_ChIP_match$nearest_gene, wildtype_ChIP_match$nearest_gene2), 'wt_CrebA_bind'] = 1
write.csv(cat_df, file = file.path(save_path, 'compiledcat_df.csv'))

###### get the scRNA logFC for late ######
save_path = file.path(TARGET_dir, 'late_data')
dir.create(save_path)
dir.create(file.path(save_path, 'plots_ct'))
dir.create(file.path(save_path, 'selected_genes_ct'))
full_gene_df = data.frame()
for(tmp_cat in unique(meta_median_2$manual_celltypes)) {
  if(tmp_cat %in% list.dirs("results/v19/DE_genes_crebA_wt", full.names =FALSE, recursive = FALSE) == FALSE) {
    next
  }
  SC_data = read.csv(file.path("results/v19/DE_genes_crebA_wt", tmp_cat, 'mut_DE_genes.csv'), row.names = 1)
  rownames(SC_data) = SC_data$feature
  sub_MA_data = MA_data[intersect(SC_data$feature, MA_data$Gene.Symbol), ]
  sub_SC_data = SC_data[intersect(SC_data$feature, MA_data$Gene.Symbol), ]
  tmp_df = data.frame('ct' = tmp_cat, 
                      'gene' = intersect(SC_data$feature, MA_data$Gene.Symbol),
                      'sc_logFC' = sub_SC_data$logFC, 
                      'ma_logFC' = sub_MA_data$logFC)
  full_gene_df = rbind(full_gene_df, tmp_df)
}

# look at the down up categories across different cell types 
interesting_genes = c()
for(tmp_ct in unique(full_gene_df$ct)) {
  subset_gene_df = full_gene_df[full_gene_df$ct == tmp_ct, ]
  subset_gene_df$categories = 'other'
  subset_gene_df[subset_gene_df$sc_logFC < -0.2 & subset_gene_df$ma_logFC < -0.3, 'categories'] = 'down'
  subset_gene_df[subset_gene_df$sc_logFC > 0.2 & subset_gene_df$ma_logFC > 0.3, 'categories'] = 'up'
  subset_gene_df[subset_gene_df$sc_logFC > -0.05 & subset_gene_df$sc_logFC < 0.05 &
                   subset_gene_df$ma_logFC > -0.05 & subset_gene_df$ma_logFC < 0.05, 'categories'] = 'static'
  p = ggplot(data = subset_gene_df, aes(x = ma_logFC, y = sc_logFC, color = categories)) +
    geom_point() +
    theme_cowplot() +
    labs(title = "Early single-cell vs MA",
         x = "MA DE Genes LogFC",
         y = "SC DE Genes logFC") + 
    ggtitle(tmp_ct)
  ggsave(file.path(save_path, 'plots_ct', paste0(tmp_ct, '.png')), plot = p, width = 7, height = 5)
  subset_gene_df = subset_gene_df[subset_gene_df$categories != 'other', ]
  interesting_genes = c(interesting_genes, subset_gene_df$gene)
  write.csv(subset_gene_df, file = file.path(save_path, 'selected_genes_ct', paste0(tmp_ct, '.csv')))
}

# get the status table across the cell types 
interesting_genes = unique(interesting_genes)
cat_df = data.frame(genes = interesting_genes)
for(tmp_ct in unique(full_gene_df$ct)) {
  subset_gene_df = full_gene_df[full_gene_df$ct == tmp_ct, ]
  subset_gene_df$categories = 'other'
  subset_gene_df[subset_gene_df$sc_logFC < -0.2 & subset_gene_df$ma_logFC < -0.3, 'categories'] = 'down'
  subset_gene_df[subset_gene_df$sc_logFC > 0.2 & subset_gene_df$ma_logFC > 0.3, 'categories'] = 'up'
  subset_gene_df[subset_gene_df$sc_logFC > -0.05 & subset_gene_df$sc_logFC < 0.05 &
                   subset_gene_df$ma_logFC > -0.05 & subset_gene_df$ma_logFC < 0.05, 'categories'] = 'static'
  rownames(subset_gene_df) = subset_gene_df$gene
  cat_df[tmp_ct] = subset_gene_df[interesting_genes, 'categories']
}
cat_df[is.na(cat_df)] = 'other'

# check if the categorized genes have CrebA binding either salivary gland specific or overall
cat_df$fkh_sage_CrebA_bind = 0
cat_df[cat_df$genes %in% c(intersect_ChIP_match$nearest_gene, intersect_ChIP_match$nearest_gene2), 'fkh_sage_CrebA_bind'] = 1
cat_df$wt_CrebA_bind = 0
cat_df[cat_df$genes %in% c(wildtype_ChIP_match$nearest_gene, wildtype_ChIP_match$nearest_gene2), 'wt_CrebA_bind'] = 1
write.csv(cat_df, file = file.path(save_path, 'compiledcat_df.csv'))