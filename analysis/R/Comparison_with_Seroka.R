library(Seurat)
library(readxl)
library(ggplot2)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "compare_with_seroka")
dir.create(TARGET_dir)

seurat_data = readRDS("accessory_data/Doe_Drosophila_Embryo_Atlas/script/curated_embryo_Doe.rds")

meta_tab = seurat_data@meta.data
meta_tab = meta_tab[meta_tab$dataset != 'stg12', ]

plot_df = data.frame(cell_types = names(table(meta_tab$cell_type)), 
                     number = as.vector(table(meta_tab$cell_type)))
plot_df$proportion = plot_df$number / sum(plot_df$number)

p<-ggplot(data=plot_df, aes(x=reorder(cell_types, proportion), y=proportion, fill = cell_types)) +
  geom_bar(stat="identity") + theme_bw() + coord_flip() + 
  ylab("Total Cell Proportion") + 
  xlab("Cell Types") + 
  ggtitle("stage 14-16: Seroka et al, 2022")
ggsave(filename = file.path(TARGET_dir, 'Seroka_cell_stg14-16_proportion.png'), plot = p, width = 14, height = 6)

# look at spearman correlation between the two 
other_data = seurat_data
our_data = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object3.rds"))

calc_spearman_correlation <- function(other_data, our_data) {
  other_meta = other_data@meta.data
  other_exp = other_data@assays$RNA@data
  
  our_meta = our_data@meta.data
  our_exp = our_data@assays$RNA@data
  
  combination_df = expand.grid(unique(our_meta$manual_celltypes), unique(other_meta$cell_type))
  colnames(combination_df) = c('our_ct', 'other_ct')  
  combination_df$spearman_correlation = NA
  
  intersect_genes = intersect(rownames(other_exp), rownames(our_exp))
  other_exp = other_exp[intersect_genes, ]
  our_exp = our_exp[intersect_genes, ]
  
  our_avg_matrix = matrix(data = NA, 
                          nrow = nrow(our_exp), 
                          ncol = length(unique(combination_df$our_ct)))
  colnames(our_avg_matrix) = unique(combination_df$our_ct)
  rownames(our_avg_matrix) = rownames(our_exp)
  
  for(ct in colnames(our_avg_matrix)) {
    sub_our_met = our_meta[our_meta$manual_celltypes == ct, ]
    sub_our_exp = our_exp[, rownames(sub_our_met)]
    our_avg = apply(sub_our_exp, FUN = mean, MARGIN = 1)
    our_avg_matrix[, ct] = our_avg
  }
  
  other_avg_matrix = matrix(data = NA, 
                          nrow = nrow(other_exp), 
                          ncol = length(unique(combination_df$other_ct)))
  colnames(other_avg_matrix) = unique(combination_df$other_ct)
  rownames(other_avg_matrix) = rownames(other_exp)
  
  for(ct in colnames(other_avg_matrix)) {
    sub_other_met = other_meta[other_meta$cell_type == ct, ]
    sub_other_exp = other_exp[, rownames(sub_other_met)]
    other_avg = apply(sub_other_exp, FUN = mean, MARGIN = 1)
    other_avg_matrix[, ct] = other_avg
  }
  
  for(temp_index in rownames(combination_df)) {
    print(temp_index)
    our_ct = combination_df[temp_index, 'our_ct']
    other_ct = combination_df[temp_index, 'other_ct']
    our_avg = our_avg_matrix[, our_ct]
    other_avg = other_avg_matrix[, other_ct]
    
    combination_df[temp_index, "spearman_correlation"] = cor(our_avg, other_avg, method = 'spearman')
  }
  
  return(combination_df)
}

library(viridis)
p = ggplot(combination_df, aes(our_ct, other_ct, fill= spearman_correlation)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) + scale_x_discrete(guide = guide_axis(angle = 45)) + 
  ggtitle("correlation similarity between our cell type and seroka's cell type")
ggsave(filename = file.path(TARGET_dir, 'seroka_cell_correlation.png'), plot = p, width = 14, height = 6)

# look at the jaccard index between top 20 marker genes 


