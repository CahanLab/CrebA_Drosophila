library(igraph)

##### load in the data ###### 
TARGET_dir = '../output/find_bound_DE_genes/'
down_DE_genes = read.csv("../output/find_bound_DE_genes/down_DE.csv", row.names = 1)
up_DE_genes = read.csv("../output/find_bound_DE_genes/up_DE.csv", row.names = 1)

# only select the down DE genes that are in sc and bound 
down_DE_genes = down_DE_genes[down_DE_genes$bound == 'True' & down_DE_genes$SC_DE == 'True', ]
up_DE_genes = up_DE_genes[up_DE_genes$bound == 'True' & up_DE_genes$SC_DE == 'True', ]

down_genes = down_DE_genes[, 'genes']
up_genes = up_DE_genes[, 'genes']

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
plot(network, vertex.color=my_color, vertex.size=3, 
     vertex.label.dist=-0.7, vertex.label.cex=0.8, 
     vertex.label.degree=-pi/2, 
     asp=0)

# Add a legend
legend("bottomleft", legend=levels(as.factor(V(network)$type))  , col = coul , bty = "n", pch=20 , pt.cex = 3, cex = 1.5, text.col=coul , horiz = FALSE, inset = c(0.1, 0.1))
