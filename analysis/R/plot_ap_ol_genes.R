TARGET_dir = file.path("results", ANALYSIS_VERSION, "Figures/ap_ol_genes")
dir.create(TARGET_dir)

###### load in the seurat objects ##### 
crebA_object = readRDS(file.path('results', ANALYSIS_VERSION, 'manual_annotation_crebA_early/manual_celltype_object.rds'))
wt_object = readRDS(file.path('results', ANALYSIS_VERSION, 'harmonized_wildtype_data/stage10-12_reharmonized_seurat.rds'))

# apodemes 
apodemes_genes = c('CG7296', 'TwdlB', 'TwdlL')

fun_color_range <- colorRampPalette(c("#f5e5e4", "#fc1303"))
my_colors <- fun_color_range(20)    

for(tmp_marker in apodemes_genes) {
  p = FeaturePlot(crebA_object, features = tmp_marker) + 
    scale_colour_gradientn(colors = my_colors)
  ggsave(filename = file.path(TARGET_dir, paste0(tmp_marker, "_crebA.png")), width = 4, height= 4)
}

# optic lobe 
ol_genes = c('E(spl)m5-HLH', 'SoxN', 'Obp99a')
fun_color_range <- colorRampPalette(c("#f5e5e4", "#fc1303"))
my_colors <- fun_color_range(20)    

for(tmp_marker in ol_genes) {
  p = FeaturePlot(wt_object, features = tmp_marker) + 
    scale_colour_gradientn(colors = my_colors)
  ggsave(filename = file.path(TARGET_dir, paste0(tmp_marker, "_wt.png")), width = 4, height= 4)
}

rownames(wt_object)[grep('m5-HLH', rownames(wt_object))]
