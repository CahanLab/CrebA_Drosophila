library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(readxl)

args = commandArgs(trailingOnly = TRUE)

##### set up the output folder path ##### 
output_path = args[3]
dir.create((output_path), recursive = TRUE)

if (grepl('4cts', args[1])[1] == TRUE) {
  color_palette = c('#84E4D9', '#AF6B7F')
  names(color_palette) = c('repression', 'activation')
} else { 
  color_palette = c('#E48482', '#6BAF92')
  names(color_palette) = c('repression', 'activation')
}

spcgs = readxl::read_excel("../input/SPCG_files/SPCG List.xlsx")
flybase_df = read.csv("../input/flybase_gene_conversion/conversion_tab.csv")
flybase_spcg_df = flybase_df[flybase_df$flybase %in% spcgs$`Drosophila FBgn`, ]

###### read in all the files ######
crebA_bound = read.csv(args[2], row.names = 1)
tss_tab = read.csv("../output/tss_table/tss_table.txt")
rownames(tss_tab) = tss_tab$gene_name

all_bound_data = data.frame()
for(tmp_index in rownames(crebA_bound)) {
  bound_loc = crebA_bound[tmp_index, 'X1'] + crebA_bound[tmp_index, 'X9']
  bound_genes = stringr::str_split(crebA_bound[tmp_index, 'nearest_gene_1'], pattern = ",")[[1]]
  bound_genes = c(bound_genes, stringr::str_split(crebA_bound[tmp_index, 'nearest_gene_2'], pattern = ",")[[1]])
  ChIP_id = crebA_bound[tmp_index, 'X3']
  bound_genes = bound_genes[bound_genes != '']
  if(length(bound_genes) == 0) {
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
include_MA = args[4]

down_DE = read.csv(file.path(args[1], "down_DE.csv"), row.names = 1)
down_DE = down_DE[down_DE$bound == 'True', ]

if(is.na(include_MA) == TRUE) {
  down_DE = down_DE[down_DE$SC_DE == 'True' | down_DE$in_situ_DE == 'True', ]
} else {
  if(include_MA == 'MA_SPCGs') {
    down_DE = down_DE[down_DE$SC_DE == 'True' | down_DE$in_situ_DE == 'True' | (down_DE$MA_DE == 'True' & down_DE$genes %in% flybase_spcg_df$gene_names), ]
  }
  else if(include_MA == 'MA_total') {
    up_DE = read.csv(file.path(args[1], "up_DE.csv"), row.names = 1)
    i_genes = intersect(down_DE$genes, up_DE$genes)
    down_DE = down_DE[down_DE$genes %in% i_genes == FALSE, ]
  } else if(include_MA == "MA_SG") {
    SG_genes = read.csv("../../analysis/results/v19/early_wt_gsea/Salivary Gland/markers_genes.csv", row.names = 1)
    SG_genes = SG_genes[SG_genes$pct.1 >= 0.1, ]
    down_DE = down_DE[down_DE$bound == 'True' & (down_DE$SC_DE == 'True' | down_DE$in_situ_DE == 'True' | (down_DE$MA_DE == 'True' & down_DE$genes %in% rownames(SG_genes))), ]
  }
}

rownames(down_DE) = down_DE$genes
sub_bound_data = all_bound_data[all_bound_data$bound_gene %in% down_DE$genes, ]

p1 = ggplot(sub_bound_data, aes(x = dist_tss)) +
  geom_histogram(aes(y = ..count../sum(..count..)), binwidth = 100, fill = color_palette[['activation']]) + 
  xlim(c(-1500, 1500)) + 
  ylab("Percentage") +
  xlab("Distance away most 5' TSS") +
  ggtitle('CrebA activated and functional genes') +
  cowplot::theme_cowplot()
ggsave(filename = file.path(output_path, 'activated_genes.png'), plot = p1, width = 5, height = 4)

##### plot out the up regulated genes ##### 
up_DE = read.csv(file.path(args[1], "up_DE.csv"), row.names = 1)
up_DE = up_DE[up_DE$bound == 'True', ]

if(is.na(include_MA) == TRUE) {
  up_DE = up_DE[up_DE$SC_DE == 'True', ]
} else {
  if(include_MA == 'MA_SPCGs') {
    up_DE = up_DE[up_DE$SC_DE == 'True' | (up_DE$MA_DE == 'True' & up_DE$genes %in% flybase_spcg_df$gene_names), ]
  }
  else if(include_MA == 'MA_total') {
    down_DE = read.csv(file.path(args[1], "down_DE.csv"), row.names = 1)
    i_genes = intersect(up_DE$genes, down_DE$genes)  
    up_DE = up_DE[up_DE$genes %in% i_genes == FALSE, ]
  }
  else if(include_MA == 'MA_SG') {
    up_DE = up_DE[up_DE$SC_DE == 'True', ]
  }
}


sub_bound_data = all_bound_data[all_bound_data$bound_gene %in% up_DE$genes, ]

p2 = ggplot(sub_bound_data, aes(x = dist_tss)) +
  geom_histogram(aes(y = ..count../sum(..count..)), binwidth = 100, fill = color_palette[['repression']]) + 
  xlim(c(-1500, 1500)) + 
  ylab("Percentage") +
  xlab("Distance from most 5' TSS") +
  ggtitle('CrebA repressed and functional genes') +
  cowplot::theme_cowplot()
ggsave(filename = file.path(output_path, 'repressed_genes.png'), plot = p2, width = 5, height = 4)

##### look at SPCGs #####
SPCG_df = readxl::read_excel("../input/SPCG_files/SPCG List.xlsx")
flybase_tab = read.csv("../input/flybase_gene_conversion/gene_names_conversion_tab.csv")
rownames(flybase_tab) = flybase_tab$flybase
sub_flybase_tab = flybase_tab[SPCG_df$`Drosophila FBgn`, ]
sub_bound_data = all_bound_data[all_bound_data$bound_gene %in% sub_flybase_tab$gene_names, ]

p3 = ggplot(sub_bound_data, aes(x = dist_tss)) +
  geom_histogram(aes(y = ..count../sum(..count..)), binwidth = 100, fill = '#87ceeb') + 
  xlim(c(-1500, 1500)) + 
  ylab("Percentage") +
  xlab("Distance from 5' TSS") +
  ggtitle('CrebA bound SPCGs') +
  cowplot::theme_cowplot()
ggsave(filename = file.path(output_path, 'bound_SPCGs_genes.png'), plot = p3, width = 5, height = 4)

##### look at how many peaks are associated with one or two ######
num_single = 0 
num_none = 0 
num_double = 0 

for (tmp_index in rownames(crebA_bound)) {
  genes = c(crebA_bound[tmp_index, 'nearest_gene_2'], crebA_bound[tmp_index, 'nearest_gene_1'])
  genes = genes[genes != '']
  if (length(genes) == 0) {
    num_none = num_none + 1
  } else if (length(genes) == 1) {
    num_single = num_single + 1
  } else if (length(genes) == 2) {
    num_double = num_double + 1
  }
}

plot_df = data.frame(ChIP_peaks = c('0', '1', '2'), 
                     number = c(num_none, num_single, num_double))

p4 = ggplot(plot_df, aes(x = ChIP_peaks, y = number)) + 
  geom_bar(stat = "identity", fill = "skyblue") +
  ylab("# of CrebA peaks") + 
  xlab("# of bound genes") + 
  theme_cowplot()

ggsave(filename = file.path(output_path, 'proportion_peaks.png'), plot = p4, width = 4, height = 4)







