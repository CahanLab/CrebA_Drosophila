TARGET_dir = file.path("results", ANALYSIS_VERSION, "Figures/plot_check_genes")
dir.create(TARGET_dir)

# load in wildtypes early data 
wt_object = readRDS("accessory_data/wild_type_seurats/GSE234602_stage10_12_seurat_object.rds")
p = FeaturePlot(wt_object, features = c('Atf6', 'CrebA', 'Met', 'Clk', 'Max'))
ggsave(filename = file.path(TARGET_dir, 'TF_genes_exp_1.png'))

p = FeaturePlot(wt_object, features = c('nau', 'Fer3', 'cato', 'crp', 'Fer1', 'sage'))
ggsave(filename = file.path(TARGET_dir, 'TF_genes_exp_2.png'))
