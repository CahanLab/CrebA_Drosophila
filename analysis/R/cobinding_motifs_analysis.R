TARGET_dir = file.path("results", ANALYSIS_VERSION, "cobinding_motifs_analysis")
dir.create(TARGET_dir)

##### get the flybase conversion #####
gtf_last_field = unique(read.table("../quantification/reference_genome_info/dmel-all-r6.33.gtf.gz", sep = "\t")[["V9"]])
gene_metadata = data.frame(
  flybase_id =
    gtf_last_field %>%
    strsplit(";") %>%
    sapply(extract2, 1) %>%
    gsub("^.*gene_id ", "", .),
  symbol =
    gtf_last_field %>%
    strsplit(";") %>%
    sapply(extract2, 2) %>%
    gsub("^.*gene_symbol ", "", .)
)
gene_metadata = gene_metadata %>% dplyr::distinct()
gene_converter = as.vector(gene_metadata$symbol)
names(gene_converter) = as.vector(gene_metadata$flybase_id)

##### get the data that ######
intersect_bind = read.csv('accessory_data/CrebA_ChIP_Analysis/output/compile_intersect_peaks_comotif/compiled_ChIP_data.csv')
intersect_bind = intersect_bind[, grep("bind", colnames(intersect_bind))]
sort(apply(intersect_bind, MARGIN = 2, FUN = sum))

wt_bind = read.csv("accessory_data/CrebA_ChIP_Analysis/output/compile_wt_peaks_comotif/compiled_ChIP_data.csv")
wt_bind = wt_bind[, grep("bind", colnames(wt_bind))]
sort(apply(wt_bind, MARGIN = 2, FUN = sum))

universe_df = data.frame(TFs_binding = colnames(intersect_bind), 
                         intersect_peaks_hits = apply(intersect_bind, MARGIN = 2, FUN = sum), 
                         wt_peaks_hits = apply(wt_bind, MARGIN = 2, FUN = sum))
write.csv(universe_df, file = file.path(TARGET_dir, 'universe_df.csv'))

universe_df$intersect_prop = universe_df$intersect_peaks_hits / nrow(intersect_bind)
universe_df$wt_peaks_prop = universe_df$wt_peaks_hits / nrow(wt_bind)
universe_df$flybase_tf = stringr::str_split_fixed(universe_df$TFs_binding, "_", n = 2)[, 1]
rownames(universe_df) = universe_df$flybase_tf
universe_df$gene_id = universe_df$flybase_tf
universe_df[rownames(universe_df), 'gene_id'] = gene_converter[rownames(universe_df)]
universe_df['FBgn0005630', 'gene_id'] = 'lola'
universe_df['FBgn0000413', 'gene_id'] = 'da'
universe_df['FBgn0000210', 'gene_id'] = 'br'
universe_df['FBgn0259172', 'gene_id'] = 'rm'

universe_df = universe_df[order(universe_df$intersect_prop, decreasing = TRUE), ]

p = ggplot(universe_df[1:10, ], aes(x = reorder(gene_id, wt_peaks_prop), y = wt_peaks_prop)) +
  geom_bar(stat="identity") + 
  theme_cowplot() + coord_flip() + 
  xlab("Co-binding Transcription Factors") +
  ylab('Proportion of CrebA Peaks with Co-binding') + 
  ggtitle("Proportion of CrebA peaks with presence of other TF motifs")
ggsave(filename = file.path(TARGET_dir, 'proportion_co_binding.png'), plot = p, height = 5, width = 5)

##### quick exploration of salivary gland cells for early #####
SG_motif_matrix = read.csv("accessory_data/CrebA_ChIP_Analysis/output/compile_intersect_peaks_comotif/co_motif_matrix.csv", row.names = 1)
bind_df = read.csv(file.path('results', ANALYSIS_VERSION, 'categorize_DE_genes/early_data/compiledcat_df.csv'), row.names = 1)
bind_df = bind_df[bind_df$Salivary.Gland != 'other', ]
bind_df = bind_df[bind_df$fkh_sage_CrebA_bind != 0, ]
pval_list = list()
for(tmp_TF in colnames(SG_motif_matrix)) {
  SG_motif_matrix = SG_motif_matrix[bind_df$genes, ]
  select_table = data.frame(down_reg = (bind_df$Salivary.Gland == 'down') * 1, 
                            TF_bind = (SG_motif_matrix[, tmp_TF] > 0) * 1)
  if(all(dim(cont_tab) == c(2, 2)) == FALSE) {
    next
  }
  cont_tab = table(select_table)
  fisher_results = fisher.test(cont_tab)
  pval_list[tmp_TF] = fisher_results$p.value
}
my_df <- t(data.frame(pval_list))

