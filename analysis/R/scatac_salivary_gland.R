library(Seurat)
library(monocle3)
library(presto)
library(pheatmap)
library(ggplot2)
set.seed(123)
TARGET_dir = file.path("results", ANALYSIS_VERSION, "scatac_salivary_gland")
dir.create(TARGET_dir)

atac_meta = readRDS("accessory_data/continuum_drosophila_embryonic_development_ATAC/output/10-16_atac_meta.rds")
atac_activity = readRDS("accessory_data/continuum_drosophila_embryonic_development_ATAC/output/10-16_activity_count_matrix.rds")

atac_meta = atac_meta[atac_meta$refined_annotation == 'Salivary gland', ]
#atac_meta = atac_meta[atac_meta$NNv1_time.new != '12-14', ]
atac_meta = atac_meta[order(atac_meta$NNv1_age), ]
select_index = c(seq(1, 90), seq(nrow(atac_meta) - 90, nrow(atac_meta)))
atac_meta = atac_meta[select_index, ]
atac_activity = atac_activity[, intersect(colnames(atac_activity), rownames(atac_meta))]
atac_meta = atac_meta[intersect(colnames(atac_activity), rownames(atac_meta)), ]
seurat_obj = Seurat::CreateSeuratObject(atac_activity)
seurat_obj@meta.data = cbind(seurat_obj@meta.data, atac_meta)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = 'LogNormalize',scale.factor = median(seurat_obj$nCount_RNA))

atac_activity_norm = seurat_obj@assays$RNA@data

results = presto::wilcoxauc(atac_activity_norm, y = atac_meta$NNv1_time.new)
results = results[results$group == '14-16', ]
rownames(results) = results$feature
#results = results[results$pval < 0.05, ]
RNA_results = read.csv(file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_salivary_gland", 'raw_startvsendtest.csv'), row.names = 1)
RNA_results = RNA_results[!is.na(RNA_results$pvalue), ]
RNA_results$adj_p = p.adjust(RNA_results$pvalue, method = 'fdr')
RNA_results$gene = rownames(RNA_results)
RNA_results = RNA_results[RNA_results$adj_p < 0.05, ]

igenes = intersect(rownames(RNA_results), rownames(results))
RNA_results = RNA_results[igenes, ]
results = results[igenes, ]

plot_df = data.frame(gene = igenes,
                     RNA_logfc = RNA_results$logFClineage1, 
                     ATAC_logfc = results$logFC)

p = ggplot(plot_df, aes(x=RNA_logfc, y=ATAC_logfc)) + geom_point() + 
  theme_bw() + ggtitle("Late vs Early diff expression and accessibility") + 
  xlab("significant RNA logFC") + 
  ylab("ATAC logFC")
ggsave(filename = file.path(TARGET_dir, "RNA_vs_ATAC_logFC.png"), plot = p, width = 10, height = 8)

TF_atac_motif = readRDS("accessory_data/continuum_drosophila_embryonic_development_ATAC/output/TF_motif_activity_matrix.rds")
TF_atac_conversion = readRDS("accessory_data/continuum_drosophila_embryonic_development_ATAC/output/TF_motif_conversion.rds")

RNA_results_TF = RNA_results[RNA_results$gene %in% TF_atac_conversion$TF_Name, ]
TF_atac_conversion = TF_atac_conversion[TF_atac_conversion$TF_Name %in% RNA_results_TF$gene, ]
TF_atac_motif = TF_atac_motif[rownames(TF_atac_motif) %in% TF_atac_conversion$Motif_ID, rownames(atac_meta)]
results = presto::wilcoxauc(TF_atac_motif, y = atac_meta$NNv1_time.new)
results = results[results$group == '14-16', ]
rownames(results) = results$feature

plot_df = data.frame()
for(region in rownames(results)) {
  TF = TF_atac_conversion[TF_atac_conversion$Motif_ID == region, 'TF_Name']
  temp_df = data.frame(region = region, 
                       TF = TF, 
                       RNA_log2fc = RNA_results_TF[TF, 'logFClineage1'], 
                       ATAC_log2fc = results[region, 'logFC'])
  plot_df = rbind(plot_df, temp_df)
}

p = ggplot(plot_df, aes(x=RNA_log2fc, y=ATAC_log2fc, color = TF)) + geom_point() + 
  theme_bw() + ggtitle("Late vs Early diff expression and TF motif score") + 
  xlab("significant RNA logFC") + 
  ylab("TF motif score logFC")

ggsave(filename = file.path(TARGET_dir, "RNA_vs_motif_logFC.png"), plot = p, width = 10, height = 8)
