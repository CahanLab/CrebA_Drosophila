library(ggplot2)
library(xml2)
library(dplyr)
library(cowplot)
library(stringr)
library(ggpubr)

###### check if the MA/SC genes are more highly expressed in embryo ###### 
dir.create("../output/reviewer_comments/embryo_overlap_correlation", recursive = TRUE)

early_obj = readRDS("../input/wild_type_seurats/GSE234602_stage10_12_seurat_object.rds")
late_obj = readRDS("../input/wild_type_seurats/GSE234602_stage13_16_seurat_object.rds")

avg_early_embryo = apply(early_obj@assays$RNA@data, MARGIN = 1, FUN = mean)
avg_late_embryo = apply(late_obj@assays$RNA@data, MARGIN = 1, FUN = mean)

down_genes = read.csv("../output/find_bound_SG_DE_genes/down_DE.csv", row.names = 1)
down_genes['early_embryo_exp'] = avg_early_embryo[down_genes$genes]
down_genes['late_embryo_exp'] = avg_late_embryo[down_genes$genes]
down_genes = down_genes[is.na(down_genes$early_embryo_exp) == FALSE, ]

down_genes['type'] = NA
down_genes[down_genes$MA_DE == 'True' & down_genes$SC_DE == 'True', 'type'] = 'MA and SC'
down_genes[down_genes$MA_DE == 'True' & down_genes$SC_DE != 'True', 'type'] = 'MA'
down_genes[down_genes$MA_DE != 'True' & down_genes$SC_DE == 'True', 'type'] = 'SC'

down_genes = down_genes[is.na(down_genes$type) == FALSE, ]

p = ggplot(down_genes, aes(x = type, y = early_embryo_exp)) + 
  geom_boxplot() + 
  theme_cowplot() + 
  xlab('CrebA activated genes') + 
  ylab("Early embryo expression") + 
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("MA", "SC"), c("MA", "MA and SC"), c('SC', 'MA and SC')), # replace with your x-axis group names
                     label = "p.format")  # shows ***, **, etc.
ggsave("../output/reviewer_comments/embryo_overlap_correlation/early_embryo.png", width = 5, height = 4)

p = ggplot(down_genes, aes(x = type, y = late_embryo_exp)) + 
  geom_boxplot() + 
  theme_cowplot() +   
  xlab('CrebA activated genes') + 
  ylab("Late embryo expression") + 
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("MA", "SC"), c("MA", "MA and SC"), c('SC', 'MA and SC')), # replace with your x-axis group names
                     label = "p.format")  # shows ***, **, etc.
ggsave("../output/reviewer_comments/embryo_overlap_correlation/late_embryo.png", width = 5, height = 4)

### plot out the boxplot without the outliers #####

p = ggplot(down_genes, aes(x = type, y = early_embryo_exp)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_cowplot() + 
  xlab('CrebA activated genes') + 
  ylab("Early embryo expression") + 
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("MA", "SC"), c("MA", "MA and SC"), c('SC', 'MA and SC')), # replace with your x-axis group names
                     label = "p.format")  # shows ***, **, etc.
ggsave("../output/reviewer_comments/embryo_overlap_correlation/early_embryo_nooutlier.png", width = 5, height = 4)

p = ggplot(down_genes, aes(x = type, y = late_embryo_exp)) + 
  geom_boxplot(outlier.shape = NA) + 
  theme_cowplot() +   
  xlab('CrebA activated genes') + 
  ylab("Late embryo expression") + 
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("MA", "SC"), c("MA", "MA and SC"), c('SC', 'MA and SC')), # replace with your x-axis group names
                     label = "p.format")  # shows ***, **, etc.
ggsave("../output/reviewer_comments/embryo_overlap_correlation/late_embryo_nooutlier.png", width = 5, height = 4)

