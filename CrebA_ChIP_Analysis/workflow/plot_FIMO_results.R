library(ggplot2)
library(xml2)
library(dplyr)
library(cowplot)
library(stringr)

##### load in required data ##### 
TF_exp_df = read.csv("../input/TF_expression_ct/all_CT_TF_percent.csv", row.names = 1)
TF_exp_df = TF_exp_df[TF_exp_df$num_cts > 20, ]

sg_TF_df = read.csv("../input/TF_expression_ct/sg_TFs.csv")
sg_TF_df = sg_TF_df[sg_TF_df$pct_in > 10, ]

sg_TF = sg_TF_df$feature
ct_TF = rownames(TF_exp_df)

# conversion dataframe 
convert_df = read.csv("../input/flybase_gene_conversion/conversion_tab.csv")
rownames(convert_df) = convert_df$flybase

sg_TF_flybase = convert_df[convert_df$gene_names %in% sg_TF, 'flybase']

##### process down regulated genes ######
input_path = file.path("../output/get_functional_peak_regions/fimo_down_DE_peaks/")
output_path = file.path(input_path, 'barplots')
color_fill = '#66c2a5'

fimo_results = read.csv(file.path(input_path, "fimo.tsv"), sep = '\t')
fimo_counts = fimo_results |> 
  group_by(motif_id) |> 
  summarise(count = n()) |> 
  arrange(desc(count))
fimo_counts$flybase = stringr::str_split_fixed(fimo_counts$motif_id, "_", n = 2)[, 1]
fimo_counts$gene_symbol = convert_df[fimo_counts$flybase, 'gene_names']
fimo_counts = fimo_counts[fimo_counts$flybase %in% sg_TF_flybase, ]

fimo_counts = fimo_counts |> 
  group_by(gene_symbol) |> 
  slice_max(count, n = 1) |> 
  arrange(desc(count))

p = ggplot(fimo_counts[1:10, ], aes(x = reorder(gene_symbol, count), y = count)) + 
  geom_bar(stat = "identity", fill = color_fill) + 
  cowplot::theme_cowplot() + 
  ylab("# Number of occurences / 117 peaks") +
  xlab("SG expressing TFs") +
  coord_flip() + 
  ggtitle(paste0('Occurences of TF motifs in activated genes CrebA peaks'))
ggsave(filename = file.path(output_path, paste0('activated_genes', ".png")), plot = p, width = 8, height = 8)

##### process down regulated genes ######
input_path = file.path("../output/get_nonfunctional_peak_regions/fimo_150_top_nonfunctional_peaks/")
output_path = file.path(input_path, 'barplots')
color_fill = '#e78ac3'

fimo_results = read.csv(file.path(input_path, "fimo.tsv"), sep = '\t')
fimo_counts = fimo_results |> 
  group_by(motif_id) |> 
  summarise(count = n()) |> 
  arrange(desc(count))
fimo_counts$flybase = stringr::str_split_fixed(fimo_counts$motif_id, "_", n = 2)[, 1]
fimo_counts$gene_symbol = convert_df[fimo_counts$flybase, 'gene_names']
fimo_counts = fimo_counts[fimo_counts$flybase %in% sg_TF_flybase, ]

fimo_counts = fimo_counts |> 
  group_by(gene_symbol) |> 
  slice_max(count, n = 1) |> 
  arrange(desc(count))

p = ggplot(fimo_counts[1:10, ], aes(x = reorder(gene_symbol, count), y = count)) + 
  geom_bar(stat = "identity", fill = color_fill) + 
  cowplot::theme_cowplot() + 
  ylab("# Number of occurences / 157 peaks") +
  xlab("SG expressing TFs") +
  coord_flip() + 
  ggtitle(paste0('Occurences of TF motifs in nonfunctional genes CrebA peaks'))
ggsave(filename = file.path(output_path, paste0('inactive_genes', ".png")), plot = p, width = 8, height = 8)

