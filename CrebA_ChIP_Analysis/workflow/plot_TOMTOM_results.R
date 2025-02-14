library(ggplot2)
library(xml2)
library(dplyr)
library(cowplot)
library(stringr)

##### load in required data ##### 

sg_TF_df = read.csv("../input/TF_expression_ct/sg_TFs.csv")
sg_TF_df = sg_TF_df[sg_TF_df$pct_in > 10, ]
sg_TF = sg_TF_df$feature

sg_DE_TF_df = sg_TF_df[sg_TF_df$logFC > 0, ]  
sg_DE_TF_df = sg_DE_TF_df[sg_DE_TF_df$pval < 0.05, ]  
sg_DE_TF = sg_DE_TF_df$feature

# conversion dataframe 
convert_df = read.csv("../input/flybase_gene_conversion/conversion_tab.csv")
rownames(convert_df) = convert_df$flybase

sg_TF_flybase = convert_df[convert_df$gene_names %in% sg_TF, 'flybase']

##### curate plotting df ###### 
curate_plot_df <- function(tomtom_results) { 
  output_list = list()
  for(tmp_seq in unique(tomtom_results$Query_ID)) {
    sub_tomtom = tomtom_results[tomtom_results$Query_ID == tmp_seq, ]
    sub_tomtom$flybase = stringr::str_split_fixed(sub_tomtom$Target_ID, "_", n = 2)[, 1]
    sub_tomtom$gene_symbol = convert_df[sub_tomtom$flybase, 'gene_names']
    sub_tomtom$logpval = -log10(sub_tomtom$p.value)
    sub_tomtom = sub_tomtom[sub_tomtom$flybase %in% sg_TF_flybase, ]
    
    sub_tomtom = sub_tomtom |> 
      group_by(flybase) |> 
      slice_max(logpval, n = 1) |> 
      ungroup() |> 
      arrange(desc(logpval))  
    sub_tomtom[, 'genes'] = paste0(sub_tomtom$flybase, "\n", sub_tomtom$gene_symbol)
    sub_tomtom[, 'symbol'] = sub_tomtom$gene_symbol
    output_list[[tmp_seq]] = sub_tomtom
  }
  return(output_list)  
}

##### process down regulated genes ######
input_path = file.path("../output/get_functional_peak_regions_SG/tomtom_down_DE_peaks/")
output_path = file.path(input_path, 'barplots')
color_fill = c('#004d00', '#66c2a5')
names(color_fill) = c('overexpressed in SG', 'expressed in SG')
tomtom_results = read.csv(file.path(input_path, "tomtom.tsv"), sep = '\t')
tomtom_results = tomtom_results[!is.na(tomtom_results$p.value), ]
plot_df_list = curate_plot_df(tomtom_results)

for(tmp_seq in names(plot_df_list)) {
  print(tmp_seq)
  plot_df = plot_df_list[[tmp_seq]]
  plot_df = plot_df[1:5, ]
  plot_df = plot_df[is.na(plot_df$Query_ID) == FALSE, ]
  plot_df$type = 'expressed in SG'
  plot_df[plot_df$symbol %in% sg_DE_TF, 'type'] = 'overexpressed in SG'
  
  p = ggplot(plot_df, aes(x = reorder(genes, logpval), y = logpval, fill = type)) + 
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = color_fill, name = 'TF type') +
    cowplot::theme_cowplot() + 
    theme(legend.position = "bottom") +
    ylab("-log10(p-val)") +
    xlab("SG expressing TFs") +
    coord_flip() + 
    ggtitle(paste0(tmp_seq, " known motifs"))
  ggsave(filename = file.path(output_path, paste0(tmp_seq, ".png")), plot = p, width = 4.5, height = 3)
}

tmp_seq = 'KACGTGK'
plot_df = plot_df_list[[tmp_seq]]
plot_df = plot_df[1:5, ]
plot_df = plot_df[is.na(plot_df$Query_ID) == FALSE, ]
plot_df$type = 'expressed in SG'
plot_df[plot_df$symbol %in% sg_DE_TF, 'type'] = 'overexpressed in SG'
p = ggplot(plot_df, aes(x = reorder(genes, logpval), y = logpval, fill = type)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = color_fill, name = 'TF type') +
  cowplot::theme_cowplot() + 
  theme(legend.position = "bottom") +
  ylab("-log10(p-val)") +
  xlab("SG expressing TFs") +
  coord_flip() + 
  ggtitle(paste0(tmp_seq, " known motifs"))
ggsave(filename = file.path(output_path, paste0(tmp_seq, "_get_template.png")), plot = p, width = 8, height = 3)
