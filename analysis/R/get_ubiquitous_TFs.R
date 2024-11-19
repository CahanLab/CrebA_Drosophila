TARGET_dir = file.path("results", ANALYSIS_VERSION, "ubiquitous_TFs")
dir.create(TARGET_dir)

flybase_TFs = read.csv("accessory_data/Drosophila_TFs/FlyBase_TFs_20241101.txt", sep = '\t', row.names = 1)
flybase_TFs$FLYBASE = rownames(flybase_TFs)

wt_obj = readRDS(file.path('results', ANALYSIS_VERSION, 'harmonized_wildtype_data/stage10-12_reharmonized_seurat.rds'))

# the minimal percent expression 
min_percent = 10
presto_output = presto::wilcoxauc(wt_obj, group_by = 'new_celltypes')
filtered_presto = presto_output[presto_output$pct_in >= min_percent, ]

filtered_presto = filtered_presto[filtered_presto$feature %in% flybase_TFs$SYMBOL, ]
filtered_presto = filtered_presto[filtered_presto$group != 'Unknown', ]
sg_TFs = filtered_presto[filtered_presto$group == 'Salivary Gland', ]

write.csv(sg_TFs, file = file.path(TARGET_dir, 'sg_TFs.csv'))

# general TF dataframe 
general_TF = matrix(data = FALSE, nrow = length(flybase_TFs$SYMBOL), ncol = length(unique(wt_obj@meta.data$new_celltypes)))
rownames(general_TF) = flybase_TFs$SYMBOL
colnames(general_TF) = unique(wt_obj@meta.data$new_celltypes)
general_TF = as.data.frame(general_TF)

for(tmp_index in rownames(filtered_presto)) {
  ct = filtered_presto[tmp_index, 'group']
  gene = filtered_presto[tmp_index, 'feature']
  general_TF[gene, ct] = TRUE
}

general_TF = general_TF[general_TF |> apply(MARGIN = 1, FUN = sum) > 0, ]
general_TF$num_cts = apply(general_TF, MARGIN = 1, FUN = sum)
general_TF = general_TF[order(general_TF$num_cts, decreasing = TRUE), ]

write.csv(general_TF, file = file.path(TARGET_dir, 'all_CT_TF_percent.csv'))
