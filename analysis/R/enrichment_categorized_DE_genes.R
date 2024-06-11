TARGET_dir = file.path("results", ANALYSIS_VERSION, "enrichment_categorized_DE_genes")
dir.create(TARGET_dir, recursive = TRUE)

categorized_DE = read.csv("accessory_data/CrebA_ChIP_Analysis/output/coordinates_bound_regions/all_bounded_genes.csv", sep = '\t')

down_genes = categorized_DE[categorized_DE$General.Type == 'Down', 'target_genes']
enriched <- enrichR::enrichr(down_genes, 'GO_Biological_Process_2018')
down_enriched_df = enriched$GO_Biological_Process_2018
write.csv(down_enriched_df, file = file.path(TARGET_dir, 'down_genes_enrichment.csv'))

target_genes = categorized_DE[categorized_DE$General.Type == 'Up', 'target_genes']
enriched <- enrichR::enrichr(target_genes, 'GO_Biological_Process_2018')
up_enriched_df = enriched$GO_Biological_Process_2018
write.csv(up_enriched_df, file = file.path(TARGET_dir, 'up_genes_enrichment.csv'))

target_genes = categorized_DE[categorized_DE$General.Type == 'Static', 'target_genes']
enriched <- enrichR::enrichr(target_genes, 'GO_Biological_Process_2018')
static_enriched_df = enriched$GO_Biological_Process_2018
write.csv(static_enriched_df, file = file.path(TARGET_dir, 'static_genes_enrichment.csv'))

down_enriched_df$logpval = -log(down_enriched_df$Adjusted.P.value)
down_enriched_df$type = 'Down'
up_enriched_df$logpval = -log(up_enriched_df$Adjusted.P.value)
up_enriched_df$type = 'Up'
static_enriched_df$logpval = -log(static_enriched_df$Adjusted.P.value)
static_enriched_df$type = 'Static'

down_enriched_df = down_enriched_df[order(down_enriched_df$Adjusted.P.value), ]
up_enriched_df = up_enriched_df[order(up_enriched_df$Adjusted.P.value), ]
static_enriched_df = static_enriched_df[order(static_enriched_df$Adjusted.P.value), ]
big_df = rbind(down_enriched_df[1:5, ], static_enriched_df[1:5, ], up_enriched_df[1:5, ])

p <- ggplot(data = big_df, aes(y = reorder(Term, logpval), x = logpval, fill = type)) +
  geom_bar(stat="identity") +
  labs(
    x = '-log10 adjusted p-value',
    y = ''
  ) + 
  scale_fill_brewer(palette = 'Set1') + 
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
    text = element_text(size = 25), 
    legend.position="none", 
    plot.title.position = "plot"
  ) + 
  theme(strip.text.x = element_blank(), axis.text.x=element_text(angle=0, vjust = 1, hjust=1)) +
  ggtitle('Geneset enrichment in down, static and up genes in CrebA mutant')  

ggsave(filename = file.path(TARGET_dir, 'enrichment_analysis.png'), plot = p, height = 10, width = 18)
