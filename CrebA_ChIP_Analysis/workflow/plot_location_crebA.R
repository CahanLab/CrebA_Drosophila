library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(readxl)

##### set up the output folder path ##### 
output_path = '../output/plot_location_crebA_binding'
dir.create((output_path), recursive = TRUE)

###### read in all the files ######
crebA_bound = read.csv('../output/match_nearest_gene/fkh_sage_intersect_genes.csv', row.names = 1)
tss_tab = read.csv("../output/tss_table/regulatory_regions_table.txt")
DE_genes = read.csv("../../analysis/results/v19/DE_genes_early_crebA_wt/Salivary Gland/mut_DE_genes.csv", row.names = 1)
rownames(DE_genes) = DE_genes$feature
rownames(tss_tab) = tss_tab$gene_name
reg_region = 1500

all_bound_data = data.frame()
for(tmp_index in rownames(crebA_bound)) {
  bound_loc = crebA_bound[tmp_index, 'X2'] - reg_region
  bound_genes = stringr::str_split(crebA_bound[tmp_index, 'in_region_gene'], pattern = ",")[[1]]
  ChIP_id = crebA_bound[tmp_index, 'X3']
  if(bound_genes[1] == '') {
    next()
  }
  for(tmp_gene in bound_genes) {
    if(tss_tab[tmp_gene, 'strand'] == "+") {
      tmp_dist = bound_loc - tss_tab[tmp_gene, 'tss']
    } else {
      tmp_dist = tss_tab[tmp_gene, 'tss'] - bound_loc
    }
    tmp_bound_data = data.frame('ChIP_id' = ChIP_id, 
                                'bound_gene' = tmp_gene, 
                                'dist_tss' = tmp_dist)
    all_bound_data = rbind(all_bound_data, tmp_bound_data)
  }
}

og_all_bound_data = all_bound_data

all_bound_data = og_all_bound_data %>%
                  group_by(bound_gene) %>%
                  slice_min(order_by = abs(dist_tss), n = 1)
write.csv(all_bound_data, file = file.path(output_path, 'all_bound_genes_dist.csv'))

###### plot out the histogram #####
down_DE = read.csv("../output/find_bound_DE_genes/down_DE.csv", row.names = 1)
down_DE = down_DE[down_DE$bound == 'True', ]
down_DE = down_DE[down_DE$SC_DE == 'True', ]
rownames(down_DE) = down_DE$genes
sub_bound_data = all_bound_data[all_bound_data$bound_gene %in% down_DE$genes, ]
sub_bound_data$logFC = abs(DE_genes[sub_bound_data$bound_gene, 'logFC'])
sub_bound_data$close_status = '>200 bp'
sub_bound_data[abs(sub_bound_data$dist_tss) < 200, 'close_status'] = "<200 bp"

p1 = ggplot(sub_bound_data, aes(x = dist_tss)) +
  geom_histogram(aes(y = ..count../sum(..count..)), binwidth = 400, fill = '#66c2a5') + 
  xlim(c(-1500, 1500)) + 
  ylim(c(0, 1)) +
  ylab("Percentage") +
  xlab("Distance away from TSS") +
  ggtitle('CrebA activated and functional genes') +
  cowplot::theme_cowplot()
ggsave(filename = file.path(output_path, 'activated_genes.png'), plot = p1, width = 5, height = 4)

p1_scatter = ggplot(sub_bound_data, aes(x = abs(dist_tss), y = logFC)) + 
  geom_point() + 
  xlab("Absolute distance") + 
  ylab("logFC") + 
  ggtitle('CrebA activated and functional genes') + 
  cowplot::theme_cowplot()
ggsave(filename = file.path(output_path, 'activated_genes_scatter.png'), plot = p1_scatter, width = 5, height = 4)

##### plot out the up regulated genes ##### 
up_DE = read.csv("../output/find_bound_DE_genes/up_DE.csv", row.names = 1)
up_DE = up_DE[up_DE$bound == 'True', ]
up_DE = up_DE[up_DE$SC_DE == 'True', ]
sub_bound_data = all_bound_data[all_bound_data$bound_gene %in% up_DE$genes, ]
sub_bound_data$logFC = abs(DE_genes[sub_bound_data$bound_gene, 'logFC'])

p2 = ggplot(sub_bound_data, aes(x = dist_tss)) +
  geom_histogram(aes(y = ..count../sum(..count..)), binwidth = 400, fill = '#fc8d62') + 
  xlim(c(-1500, 1500)) + 
  ylim(c(0, 1)) +
  ylab("Percentage") +
  xlab("Distance away from TSS") +
  ggtitle('CrebA repressed and functional genes') +
  cowplot::theme_cowplot()
ggsave(filename = file.path(output_path, 'repressed_genes.png'), plot = p2, width = 5, height = 4)

p2_scatter = ggplot(sub_bound_data, aes(x = abs(dist_tss), y = logFC)) + 
  geom_point() + 
  xlab("Absolute distance") + 
  ylab("logFC") + 
  ggtitle('CrebA repressed and functional genes') + 
  cowplot::theme_cowplot()
ggsave(filename = file.path(output_path, 'repressed_genes_scatter.png'), plot = p2_scatter, width = 5, height = 4)

##### look at SPCGs #####
SPCG_df = readxl::read_excel("../input/SPCG_files/SPCG List.xlsx")
flybase_tab = read.csv("../input/flybase_gene_conversion/gene_names_conversion_tab.csv")
rownames(flybase_tab) = flybase_tab$flybase
sub_flybase_tab = flybase_tab[SPCG_df$`Drosophila FBgn`, ]
sub_bound_data = all_bound_data[all_bound_data$bound_gene %in% sub_flybase_tab$gene_names, ]

p3 = ggplot(sub_bound_data, aes(x = dist_tss)) +
  geom_histogram(aes(y = ..count../sum(..count..)), binwidth = 400, fill = '#87ceeb') + 
  xlim(c(-1500, 1500)) + 
  ylim(c(0, 1)) +
  ylab("Percentage") +
  xlab("Distance away from TSS") +
  ggtitle('CrebA bound SPCGs') +
  cowplot::theme_cowplot()
ggsave(filename = file.path(output_path, 'bound_SPCGs_genes.png'), plot = p3, width = 5, height = 4)
