library(enrichR)
library(ggplot2)
library(igraph)
library(readxl)
library(stringr)

enrichR::setEnrichrSite('FlyEnrichr')

args = commandArgs(trailingOnly = TRUE)

input_dir = file.path(args[1])
TARGET_dir = file.path(args[2])
include_MA = args[3]

dir.create(TARGET_dir, recursive = TRUE)

# load in spcgs 
spcgs = readxl::read_excel("../input/SPCG_files/SPCG List.xlsx")
flybase_df = read.csv("../input/flybase_gene_conversion/conversion_tab.csv")
flybase_spcg_df = flybase_df[flybase_df$flybase %in% spcgs$`Drosophila FBgn`, ]

##### get the DE genes ###### 
down_DE_genes = read.csv(file.path(input_dir, "down_DE.csv"), row.names = 1)
up_DE_genes = read.csv(file.path(input_dir, "up_DE.csv"), row.names = 1)

if (grepl('4cts', input_dir)[1] == TRUE) {
  color_palette = c('#84E4D9', '#AF6B7F')
  names(color_palette) = c('repression', 'activation')
  num_entry = 4
} else { 
  color_palette = c('#E48482', '#6BAF92')
  names(color_palette) = c('repression', 'activation')
  num_entry = 5
}

# only select the down DE genes that are in sc and bound 
if(is.na(include_MA) == TRUE) {
  down_DE_genes = down_DE_genes[down_DE_genes$bound == 'True' & (down_DE_genes$SC_DE == 'True' | down_DE_genes$in_situ_DE == 'True'), ]
  up_DE_genes = up_DE_genes[up_DE_genes$bound == 'True' & up_DE_genes$SC_DE == 'True', ]
} else {
  if(include_MA == 'MA_SPCGs') {
    down_DE_genes = down_DE_genes[down_DE_genes$bound == 'True' & (down_DE_genes$SC_DE == 'True' | down_DE_genes$in_situ_DE == 'True' | (down_DE_genes$MA_DE == 'True' & down_DE_genes$genes %in% flybase_spcg_df$gene_names)), ]
    up_DE_genes = up_DE_genes[up_DE_genes$bound == 'True' & (up_DE_genes$SC_DE == 'True' | (up_DE_genes$MA_DE == 'True' & up_DE_genes$genes %in% flybase_spcg_df$gene_names)), ]
  }
  else if(include_MA == 'MA_total') {
    i_genes = intersect(down_DE_genes$genes, up_DE_genes$genes)
    print(i_genes)
    down_DE_genes = down_DE_genes[down_DE_genes$genes %in% i_genes == FALSE, ]
    up_DE_genes = up_DE_genes[up_DE_genes$genes %in% i_genes == FALSE, ]
  } else if(include_MA == "MA_SG") {
    SG_genes = read.csv("../../analysis/results/v19/early_wt_gsea/Salivary Gland/markers_genes.csv", row.names = 1)
    SG_genes = SG_genes[SG_genes$pct.1 >= 0.1, ]
    
    # load in the diff expressed genes between mutants and wildtype
    mut_diff_genes = read.csv("../../analysis/results/v19/DE_genes_early_crebA_wt/Salivary Gland/mut_DE_genes.csv", row.names = 1)
    mut_diff_genes = mut_diff_genes[mut_diff_genes$logFC < 0, ]
    MA_select_genes = intersect(rownames(SG_genes), mut_diff_genes$feature)
    
    down_DE_genes = down_DE_genes[down_DE_genes$bound == 'True' & (down_DE_genes$SC_DE == 'True' | down_DE_genes$in_situ_DE == 'True' | (down_DE_genes$MA_DE == 'True' & down_DE_genes$genes %in% MA_select_genes)), ]
    up_DE_genes = up_DE_genes[up_DE_genes$bound == 'True' & (up_DE_genes$SC_DE == 'True'), ]
  }
}

##### categorize them ##### 
down_genes = down_DE_genes[, 'genes']
enriched <- enrichR::enrichr(down_genes, 'GO_Biological_Process_2018')
down_enriched_df = enriched$GO_Biological_Process_2018
write.csv(down_enriched_df, file = file.path(TARGET_dir, 'down_genes_enrichment.csv'))

up_genes = up_DE_genes[, 'genes']
enriched <- enrichR::enrichr(up_genes, 'GO_Biological_Process_2018')
up_enriched_df = enriched$GO_Biological_Process_2018
write.csv(up_enriched_df, file = file.path(TARGET_dir, 'up_genes_enrichment.csv'))

##### look at the up and down genes #####
down_enriched_df$logpval = -log10(down_enriched_df$Adjusted.P.value)
down_enriched_df$type = 'Activation'
up_enriched_df$logpval = -log10(up_enriched_df$Adjusted.P.value)
up_enriched_df$type = 'Repression'

down_enriched_df = down_enriched_df[order(down_enriched_df$Adjusted.P.value), ]
up_enriched_df = up_enriched_df[order(up_enriched_df$Adjusted.P.value), ]
big_df = rbind(down_enriched_df[1:num_entry, ], up_enriched_df[1:num_entry, ])

for(tmp_row in rownames(big_df)) {
  tmp_term = big_df[tmp_row, 'Term']  
  tmp_list = stringr::str_split(tmp_term, " ")[[1]]
  if(length(tmp_list) >= 10) {
    tmp_list = append(tmp_list, "\n", after = round(length(tmp_list)/2))
    new_name = paste(unlist(tmp_list), collapse = " ")
    big_df[tmp_row, 'Term'] = new_name
  }
}

p <- ggplot(data = big_df, aes(y = reorder(Term, logpval), x = logpval, fill = type)) +
  geom_bar(stat="identity") +
  labs(
    x = '-log10 adjusted p-value',
    y = ''
  ) + 
  scale_fill_manual(values = c("Repression" = color_palette[['repression']], "Activation" = color_palette[['activation']])) +
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
    text = element_text(size = 35), 
    legend.position="none", 
    plot.title.position = "plot"
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=0, vjust = 1, hjust=1)) +
  ggtitle('Geneset enrichment of CrebA bound and functional genes')  

if(include_MA == 'MA_total') {
  ggsave(filename = file.path(TARGET_dir, 'enrichment_analysis.png'), plot = p, height = 10, width = 24)

} else {
  ggsave(filename = file.path(TARGET_dir, 'enrichment_analysis.png'), plot = p, height = 10, width = 18)
}

##### remove the spcgs ##### 
down_genes = down_genes[down_genes %in% flybase_spcg_df$gene_names == FALSE]

enriched <- enrichR::enrichr(down_genes, 'GO_Biological_Process_2018')
down_enriched_df = enriched$GO_Biological_Process_2018
write.csv(down_enriched_df, file = file.path(TARGET_dir, 'down_genes_no_spcg_enrichment.csv'))

up_genes = up_DE_genes[, 'genes']
enriched <- enrichR::enrichr(up_genes, 'GO_Biological_Process_2018')
up_enriched_df = enriched$GO_Biological_Process_2018
write.csv(up_enriched_df, file = file.path(TARGET_dir, 'up_genes_no_spcg_enrichment.csv'))

##### look at the up and down genes #####
down_enriched_df$logpval = -log10(down_enriched_df$Adjusted.P.value)
down_enriched_df$type = 'Activation'
up_enriched_df$logpval = -log10(up_enriched_df$Adjusted.P.value)
up_enriched_df$type = 'Repression'

down_enriched_df = down_enriched_df[order(down_enriched_df$Adjusted.P.value), ]
up_enriched_df = up_enriched_df[order(up_enriched_df$Adjusted.P.value), ]
big_df = rbind(down_enriched_df[1:num_entry, ], up_enriched_df[1:num_entry, ])

for(tmp_row in rownames(big_df)) {
  tmp_term = big_df[tmp_row, 'Term']  
  tmp_list = stringr::str_split(tmp_term, " ")[[1]]
  if(length(tmp_list) >= 10) {
    tmp_list = append(tmp_list, "\n", after = round(length(tmp_list)/2))
    new_name = paste(unlist(tmp_list), collapse = " ")
    big_df[tmp_row, 'Term'] = new_name
  }
}

p <- ggplot(data = big_df, aes(y = reorder(Term, logpval), x = logpval, fill = type)) +
  geom_bar(stat="identity") +
  labs(
    x = '-log10 adjusted p-value',
    y = ''
  ) + 
  scale_fill_manual(values = c("Repression" = color_palette[['repression']], "Activation" = color_palette[['activation']])) +
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
    text = element_text(size = 35), 
    legend.position="none", 
    plot.title.position = "plot"
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=0, vjust = 1, hjust=1)) +
  ggtitle('Geneset enrichment of CrebA bound and functional genes (no SPCGs)')  

if(include_MA == 'MA_total') {
  ggsave(filename = file.path(TARGET_dir, 'enrichment_no_spcg_analysis.png'), plot = p, height = 10, width = 24)
} else {
  ggsave(filename = file.path(TARGET_dir, 'enrichment_no_spcg_analysis.png'), plot = p, height = 10, width = 19)
}
