library(ggplot2)
library(dplyr)
library(cowplot)

##### create the directory ######
TARGET_dir = "../output/motif_counts/"
dir.create(TARGET_dir)

ACGTG_path = '../output/ACGTG_cores/'
ACGTC_path = '../output/ACGTC_cores'

file_to_type = list()
file_to_type[['down_DE_peaks.csv']] = 'activation'
file_to_type[['up_DE_peaks.csv']] = 'repression'
file_to_type[['150_top_nonfunctional_peaks.csv']] = 'non-functional'

# this is to compile the motifs
compile_motifs <- function(file_name = 'down_DE_peaks.csv') { 
  ACGTG_df = read.csv(file.path(ACGTG_path, file_name), row.names = 1)
  colnames(ACGTG_df) = paste0("ACGTG_", colnames(ACGTG_df))
  
  ACGTC_df = read.csv(file.path(ACGTC_path, file_name), row.names = 1)
  colnames(ACGTC_df) = paste0("ACGTC_", colnames(ACGTC_df))
  
  combined_df = cbind(ACGTG_df, ACGTC_df)
  return(combined_df)
}

for(tmp_file in c('down_DE_peaks.csv', 'up_DE_peaks.csv', '150_top_nonfunctional_peaks.csv')) {
  combined_df = compile_motifs(tmp_file)
  write.csv(combined_df, file.path(TARGET_dir, paste0(file_to_type[[tmp_file]], '.csv')))
}

tmp_matrix = matrix(data = NA, nrow = 3, ncol = 3)
rownames(tmp_matrix) = c('activation', 'repression', 'non-functional')
colnames(tmp_matrix) = c('acgtg', 'acgtc', 'either')

for(tmp_type in c('activation', 'repression', 'non-functional')) {
  tmp_df = read.csv(file.path(TARGET_dir, paste0(tmp_type, '.csv')), row.names = 1)
  tmp_matrix[tmp_type, 'acgtg'] = sum(tmp_df$ACGTG_motif != "" | tmp_df$ACGTG_comp_rev_motif != "") / nrow(tmp_df)
  tmp_matrix[tmp_type, 'acgtc'] = sum(tmp_df$ACGTC_motif != "" | tmp_df$ACGTC_comp_rev_motif != "") / nrow(tmp_df)
  tmp_matrix[tmp_type, 'either'] = sum(tmp_df$ACGTG_motif != "" | tmp_df$ACGTG_comp_rev_motif != "" | tmp_df$ACGTC_motif != "" | tmp_df$ACGTC_comp_rev_motif != "") / nrow(tmp_df)
  
}

plot_df = data.frame()
for(tmp_row in rownames(tmp_matrix)) { 
  tmp_plot_df = data.frame(gene_type = tmp_row, 
                           motif_type = colnames(tmp_matrix), 
                           percentage = tmp_matrix[tmp_row, ]) 
  plot_df = rbind(plot_df, tmp_plot_df)
}

p = ggplot(plot_df, aes(x = motif_type, y = percentage, fill = gene_type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("#6BAF92", "#A4A9D1", "#E48482")) +  # Custom colors
  labs(
    title = "",
    x = "Motif cores",
    y = "Percentage of CrebA binding with motif core",
    fill = "Target gene type"
  ) +
  cowplot::theme_cowplot()
ggsave(filename = file.path(TARGET_dir, 'barplot.png'), plot = p, width = 8, height = 5)
