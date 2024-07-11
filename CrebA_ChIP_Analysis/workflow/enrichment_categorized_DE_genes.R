library(enrichR)
library(ggplot2)
library(igraph)

enrichR::setEnrichrSite('FlyEnrichr')

TARGET_dir = file.path("../output/enrichment_categorized_DE_genes")
dir.create(TARGET_dir, recursive = TRUE)

##### get the DE genes ###### 
down_DE_genes = read.csv("../output/find_bound_DE_genes/down_DE.csv", row.names = 1)
up_DE_genes = read.csv("../output/find_bound_DE_genes/up_DE.csv", row.names = 1)

# only select the down DE genes that are in sc and bound 
down_DE_genes = down_DE_genes[down_DE_genes$bound == 'True' & down_DE_genes$SC_DE == 'True', ]
up_DE_genes = up_DE_genes[up_DE_genes$bound == 'True' & up_DE_genes$SC_DE == 'True', ]

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
down_enriched_df$logpval = -log(down_enriched_df$Adjusted.P.value)
down_enriched_df$type = 'Activation'
up_enriched_df$logpval = -log(up_enriched_df$Adjusted.P.value)
up_enriched_df$type = 'Repression'

down_enriched_df = down_enriched_df[order(down_enriched_df$Adjusted.P.value), ]
up_enriched_df = up_enriched_df[order(up_enriched_df$Adjusted.P.value), ]
big_df = rbind(down_enriched_df[1:5, ], up_enriched_df[1:5, ])

p <- ggplot(data = big_df, aes(y = reorder(Term, logpval), x = logpval, fill = type)) +
  geom_bar(stat="identity") +
  labs(
    x = '-log10 adjusted p-value',
    y = ''
  ) + 
  scale_fill_brewer(palette = 'Set2') + 
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
  ggtitle('Geneset enrichment of CrebA bound and functional genes')  

ggsave(filename = file.path(TARGET_dir, 'enrichment_analysis.png'), plot = p, height = 10, width = 18)

###### plot out the network for downstream genes of CrebA ######
library(igraph)
links = data.frame(
  'source' = rep(x = 'CrebA', length(c(down_genes, up_genes))), 
  'target' = c(down_genes, up_genes)
)
nodes = data.frame(
  names = c(down_genes[down_genes != 'CrebA'], up_genes, 'CrebA'), 
  type = c(rep('Activation', length(down_genes[down_genes != 'CrebA'])), rep('Repression', length(up_genes)), 'CrebA')
)

# Turn it into igraph object
network <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 

# Make a palette of 3 colors
library(RColorBrewer)
coul  <- brewer.pal(3, "Set2") 

# Create a vector of color
V(network)$type = factor(V(network)$type, levels = c('Activation', 'Repression', 'CrebA'))
my_color <- coul[as.numeric(as.factor(V(network)$type))]
E(network)$width = 0.1

# Make the plot
png(filename = file.path(TARGET_dir, 'creba_grn.png'), width = 5, height = 5)
plot(network, vertex.color=my_color, vertex.size=3, 
     vertex.label.dist=-0.7, vertex.label.cex=0.8, 
     vertex.label.degree=-pi/2, 
     asp=0)
dev.off()

# Add a legend
legend("bottomleft", legend=levels(as.factor(V(network)$type))  , col = coul , bty = "n", pch=20 , pt.cex = 3, cex = 1.5, text.col=coul , horiz = FALSE, inset = c(0.1, 0.1))
