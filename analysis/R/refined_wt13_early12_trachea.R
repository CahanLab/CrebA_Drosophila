library(Seurat)
library(monocle3)
library(presto)
library(enrichR)

set.seed(123)
TARGET_dir = file.path("results", ANALYSIS_VERSION, "refined_wt_late_early_trachea")
dir.create(TARGET_dir)

wt_late_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_wt13/manual_celltype_object3.rds"))

##########################
# TODO the below will be changed as time goes on 
wt_early_object = readRDS(file.path("results", ANALYSIS_VERSION, "early_wt12_integrated/BDGP_automated_annotation_object.rds"))
sub_wt_early_object = subset(wt_early_object, subset = Integrated_tentativeCellType == 'tracheal primordium')
sub_wt_early_object$experimental_condition = 'early'
########################

sub_wt_late_object = subset(wt_late_object, subset = manual_celltypes == "Trachea")
sub_wt_late_object$experimental_condition = 'late'

combined_ct_object = merge(sub_wt_early_object, sub_wt_late_object)
object = Seurat::CreateSeuratObject(combined_ct_object@assays$RNA@counts, project = 'Trachea')
object@meta.data$experimental_condition = combined_ct_object@meta.data$experimental_condition
object@meta.data$batch = paste0(combined_ct_object@meta.data$experimental_condition, "_", combined_ct_object@meta.data$batch)

object = Seurat::NormalizeData(object)
cellCycleMarkers = read.csv("accessory_data/cellCycleMarkers.csv", skip = 1, header = T)
object %<>% CellCycleScoring(s.features = cellCycleMarkers$S.phase.markers., g2m.features = cellCycleMarkers$G2.M.phase.markers.)

# start monocle3 
expression_matrix = object@assays$RNA@counts
cell_metadata = object@meta.data
gene_annotation = object@assays$RNA@meta.features
gene_annotation$gene_short_name = rownames(gene_annotation)
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 100)

cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=TRUE,  color_cells_by = "batch", cell_size = 1, label_cell_groups = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_UMAP.png'), width = 8, height = 6)

plot_cells(cds,
           genes=c("Osi6"),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE, cell_size = 2)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_toe.png'), width = 8, height = 6)

cds <- cluster_cells(cds, resolution = 1e-2)
plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE, cell_size = 1)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_cluster.png'), width = 8, height = 6)

marker_test_res <- monocle3::top_markers(cds, group_cells_by="cluster", 
                               reference_cells=1000, cores=8)

write.csv(marker_test_res, file = file.path(TARGET_dir, "top_markers_monocle.csv"))
rank_sum_results = presto::wilcoxauc(normalized_counts(cds), cds@clusters$UMAP$clusters)
write.csv(rank_sum_results, file = file.path(TARGET_dir, "rank_sum_test.csv"))

plot_cells(cds, color_cells_by = "Phase", label_cell_groups = FALSE, cell_size = 1, show_trajectory_graph = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_phase.png'), width = 8, height = 6)

cds <- learn_graph(cds)

plot_cells(cds,
           color_cells_by = "batch",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
cds <- order_cells(cds)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5, cell_size = 1, show_trajectory_graph = FALSE)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_pt.png'), width = 8, height = 6)
saveRDS(cds, file = file.path(TARGET_dir, "monocle3_no_batch_correct_object.rds"))

# this is where we use TradeSeq to test association 
expRaw = expression_matrix[apply(expression_matrix, MARGIN = 1, FUN = function(x){return(sum(x > 0) > 10) }), ]
pt<-as.data.frame(pseudotime(cds))
colnames(pt)<-"pseudotime"
cw<-as.matrix(rep(1,nrow(pt)))
rownames(cw)<-rownames(pt)
ts<-tradeSeq::fitGAM(as.matrix(expRaw),pseudotime=as.matrix(pt),cellWeights=cw, parallel = TRUE)
saveRDS(ts, file = file.path(TARGET_dir, "tradeseq_fitgam_results.rds"))

ATres<-associationTest(ts)
saveRDS(ATres, file = file.path(TARGET_dir, 'raw_associationTest.rds'))

ATres = ATres[!is.na(ATres$pvalue), ]
ATres$adj_p = p.adjust(ATres$pvalue)
ATres = ATres[ATres$adj_p < 0.05, ]
write.csv(ATres, file = file.path(TARGET_dir, "significant_associationTest.csv"))

startRes <- startVsEndTest(ts)
startRes$adj_p = p.adjust(startRes$pvalue)
write.csv(startRes, file = file.path(TARGET_dir, 'raw_startvsendtest.csv'))  

##################################################
# find the identity of the weird cluster 
# figure out the weird branch 
cds = readRDS(file.path(TARGET_dir, "monocle3_no_batch_correct_object.rds"))
plot_cells(cds, color_cells_by = "cluster", 
           label_cell_groups = TRUE, 
           cell_size = 1, show_trajectory_graph = FALSE, group_label_size = 10)
ggsave(file.path(TARGET_dir, 'monocle3_no_batch_corrected_clusters.png'), width = 8, height = 6)

# get the biological process that are enriched in each cluster 
dir.create(file.path(TARGET_dir, 'cluster_process'))
marker_test_res = read.csv(file.path(TARGET_dir, "top_markers_monocle.csv"), row.names = 1)
withr::with_dir(file.path(TARGET_dir, 'cluster_process'), {
  for(cluster in unique(marker_test_res$cell_group)) {
    sub_marker_test_res = marker_test_res[marker_test_res$cell_group == cluster, ]
    enrichment_results = enrichR::enrichr(
      genes = sub_marker_test_res$gene_id, 
      databases = c(
        "GO_Biological_Process_2018", 
        "GO_Molecular_Function_2018"
      )
    )
    biological_analysis = enrichment_results$GO_Biological_Process_2018
    write.csv(biological_analysis, file = paste0(cluster, "_biological_process.csv"))
  }
})

withr::with_dir(file.path(TARGET_dir, 'cluster_process'), {
  sub_marker_test_res = marker_test_res[marker_test_res$cell_group == 10, ]
  plot_cells(cds,
             genes=sub_marker_test_res$gene_id,
             label_cell_groups=FALSE,
             show_trajectory_graph=FALSE, cell_size = 2)
  ggsave(filename = 'cluster_10_genes.png')
})


library(enrichR)
cds = readRDS(file.path(TARGET_dir, "monocle3_no_batch_correct_object.rds"))
enrichR::setEnrichrSite("FlyEnrichr")
startRes = read.csv(file.path(TARGET_dir, 'raw_startvsendtest.csv'), row.names = 1)
startRes = startRes[!is.na(startRes$pvalue), ]
startRes$adj_p = p.adjust(startRes$pvalue, method = 'fdr')
startRes = startRes[startRes$adj_p < 0.05, ]

startRes_early = startRes[startRes$logFClineage1 < 0, ]
startRes_late = startRes[startRes$logFClineage1 >= 0, ]

enrichment_results = enrichR::enrichr(
  genes = rownames(startRes_early), 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018"
  )
)

biological_analysis = enrichment_results$GO_Biological_Process_2018
molecular_analysis = enrichment_results$GO_Molecular_Function_2018
write.csv(biological_analysis, file = file.path(TARGET_dir, 'sig_GO_biological_early.csv'))

enrichment_results = enrichR::enrichr(
  genes = rownames(startRes_late), 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018"
  )
)

biological_analysis = enrichment_results$GO_Biological_Process_2018
molecular_analysis = enrichment_results$GO_Molecular_Function_2018
write.csv(biological_analysis, file = file.path(TARGET_dir, 'sig_GO_biological_late.csv'))

#########################################################################
biological_analysis = read.csv(file.path(TARGET_dir, 'sig_GO_biological_early.csv'), row.names = 1)
biological_analysis = biological_analysis[order(biological_analysis$P.value), ]
sub_biological_analysis = biological_analysis[1:15, ]
sub_biological_analysis$log_p = -log10(sub_biological_analysis$Adjusted.P.value)
p <- ggplot(data=sub_biological_analysis, aes(x=reorder(Term, log_p), y=log_p)) +
  xlab("GO Biological Process") +
  ylab("-log10 P-value") +
  ggtitle("Early Genes Enrichment") +
  geom_bar(stat="identity", fill="steelblue") + coord_flip() + theme_bw()
ggsave(filename = file.path(TARGET_dir, "early_go_enrichment.png"), plot = p, width = 10, height = 8)

biological_analysis = read.csv(file.path(TARGET_dir, 'sig_GO_biological_late.csv'), row.names = 1)
biological_analysis = biological_analysis[order(biological_analysis$Adjusted.P.value), ]
sub_biological_analysis = biological_analysis[1:15, ]
sub_biological_analysis$log_p = -log10(sub_biological_analysis$Adjusted.P.value)
p <- ggplot(data=sub_biological_analysis, aes(x=reorder(Term, log_p), y=log_p)) +
  xlab("GO Biological Process") +
  ylab("-log10 adjusted P-value") +
  ggtitle("Late Genes Enrichment") +
  geom_bar(stat="identity", fill="salmon") + coord_flip() + theme_bw()
ggsave(filename = file.path(TARGET_dir, "late_go_enrichment.png"), plot = p, width = 10, height = 8)

# early genes 
dir.create(file.path(TARGET_dir, 'early_genes'))
biological_analysis = read.csv(file.path(TARGET_dir, 'sig_GO_biological_early.csv'), row.names = 1)
biological_analysis = biological_analysis[order(biological_analysis$P.value), ]
biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
biological_analysis = biological_analysis[order(biological_analysis$Combined.Score, decreasing = TRUE), ]
biological_analysis = biological_analysis[1:50, ]
withr::with_dir(file.path(TARGET_dir, 'early_genes'), {
  for(term in biological_analysis$Term) {
    target_genes = biological_analysis[biological_analysis$Term == term, 'Genes']
    term = stringr::str_replace(term, "/", "-")
    target_genes = stringr::str_split(target_genes, ";")[[1]]
    #target_genes = stringr::str_to_title(target_genes)
    norm_exp = monocle3::normalized_counts(cds)
    cluster_id = monocle3::clusters(cds)
    norm_exp = as.matrix(norm_exp)
    norm_exp = norm_exp[target_genes, ]
    
    # this will change 
    #norm_exp = norm_exp[apply(norm_exp, MARGIN = 1, FUN = max) > 1, ]
    pt = monocle3::pseudotime(cds)
    pt = data.frame(pseudotime = pt, 
                    cluster_id = monocle3::clusters(cds))
    pt = pt[pt$cluster_id != 10, ]
    pt$cluster_id = NULL
    plot_df = cbind(pt, t(norm_exp[, rownames(pt)]))
    smoothed_df = data.frame()
    for(gene in colnames(plot_df)) {
      if(gene == 'pseudotime') {
        next
      }
      else {
        yy = ksmooth(plot_df[, 'pseudotime'], plot_df[, gene], kernel="normal", bandwidth = 1.5, x.points=plot_df[, 'pseudotime'])
        if(nrow(smoothed_df) == 0) {
          smoothed_df = data.frame('pseudotime' = yy$x)
        }
        smoothed_df[, gene] = yy$y
      }
    }
    smoothed_df$pseudotime = NULL
    smoothed_df = t(smoothed_df)
    scaled_exp = t(scale(t(smoothed_df)))
    sorted_genes = names(sort(apply(scaled_exp, MARGIN = 1, FUN = which.max)))
    scaled_exp = scaled_exp[sorted_genes, ]
    pdf(file = paste0(term, "_dynamic_gene_heatmap.pdf"), height = 10)
    pheatmap(scaled_exp[sorted_genes, ], cluster_cols = FALSE, cluster_rows = FALSE)
    dev.off()
    
    convert_line_plot <- function(scaled_exp) {
      plot_df = data.frame()
      for(gene in rownames(scaled_exp)) {
        temp_plot = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                               scaled_exp = scaled_exp[gene, ], 
                               gene = gene)
        plot_df = rbind(plot_df, temp_plot)
      }
      
      plot_df$pseudotime = plot_df$pseudotime / max(plot_df$pseudotime)
      return(plot_df)
    }
    
    plot_df = convert_line_plot(scaled_exp)
    p<-ggplot(plot_df, aes(x=pseudotime, y=scaled_exp, group=gene)) +
      xlab("pseudotime") + 
      ylab("smoothed and scaled expression") +
      ggtitle(paste0('genes in ', term)) +
      geom_line(aes(color=gene)) + theme_bw() 
    
    ggsave(paste0(term, "_dynamic_gene_line.png"), plot = p, width = 10, height = 8)
  }
})

# early genes 
dir.create(file.path(TARGET_dir, 'late_genes'))
biological_analysis = read.csv(file.path(TARGET_dir, 'sig_GO_biological_late.csv'), row.names = 1)
biological_analysis = biological_analysis[order(biological_analysis$P.value), ]
biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
biological_analysis = biological_analysis[order(biological_analysis$Combined.Score, decreasing = TRUE), ]
biological_analysis = biological_analysis[1:50, ]
withr::with_dir(file.path(TARGET_dir, 'late_genes'), {
  for(term in biological_analysis$Term) {
    target_genes = biological_analysis[biological_analysis$Term == term, 'Genes']
    term = stringr::str_replace(term, "/", "-")
    target_genes = stringr::str_split(target_genes, ";")[[1]]
    #target_genes = stringr::str_to_title(target_genes)
    norm_exp = monocle3::normalized_counts(cds)
    cluster_id = monocle3::clusters(cds)
    norm_exp = as.matrix(norm_exp)
    norm_exp = norm_exp[target_genes, ]
    
    # this will change 
    #norm_exp = norm_exp[apply(norm_exp, MARGIN = 1, FUN = max) > 1, ]
    pt = monocle3::pseudotime(cds)
    pt = data.frame(pseudotime = pt, 
                    cluster_id = monocle3::clusters(cds))
    pt = pt[pt$cluster_id != 10, ]
    pt$cluster_id = NULL
    plot_df = cbind(pt, t(norm_exp[, rownames(pt)]))
    smoothed_df = data.frame()
    for(gene in colnames(plot_df)) {
      if(gene == 'pseudotime') {
        next
      }
      else {
        yy = ksmooth(plot_df[, 'pseudotime'], plot_df[, gene], kernel="normal", bandwidth = 1.5, x.points=plot_df[, 'pseudotime'])
        if(nrow(smoothed_df) == 0) {
          smoothed_df = data.frame('pseudotime' = yy$x)
        }
        smoothed_df[, gene] = yy$y
      }
    }
    smoothed_df$pseudotime = NULL
    smoothed_df = t(smoothed_df)
    scaled_exp = t(scale(t(smoothed_df)))
    sorted_genes = names(sort(apply(scaled_exp, MARGIN = 1, FUN = which.max)))
    scaled_exp = scaled_exp[sorted_genes, ]
    pdf(file = paste0(term, "_dynamic_gene_heatmap.pdf"), height = 10)
    pheatmap(scaled_exp[sorted_genes, ], cluster_cols = FALSE, cluster_rows = FALSE)
    dev.off()
    
    convert_line_plot <- function(scaled_exp) {
      plot_df = data.frame()
      for(gene in rownames(scaled_exp)) {
        temp_plot = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                               scaled_exp = scaled_exp[gene, ], 
                               gene = gene)
        plot_df = rbind(plot_df, temp_plot)
      }
      
      plot_df$pseudotime = plot_df$pseudotime / max(plot_df$pseudotime)
      return(plot_df)
    }
    
    plot_df = convert_line_plot(scaled_exp)
    p<-ggplot(plot_df, aes(x=pseudotime, y=scaled_exp, group=gene)) +
      xlab("pseudotime") + 
      ylab("smoothed and scaled expression") +
      ggtitle(paste0('genes in ', term)) +
      geom_line(aes(color=gene)) + theme_bw() 
    
    ggsave(paste0(term, "_dynamic_gene_line.png"), plot = p, width = 10, height = 8)
  }
})

####################################################
# rib_genes
excel_path = file.path("accessory_data/ChIP_Targets/SG-and-Tracheal-Rib_IDR_20220928.xlsx")
page_names = openxlsx::getSheetNames(excel_path)
ChIP_genes = openxlsx::read.xlsx(excel_path, sheet = page_names[2], colNames = TRUE)
dir.create(file.path(TARGET_dir, 'rib_genes'))

biological_analysis = read.csv(file.path(TARGET_dir, 'sig_GO_biological_late.csv'), row.names = 1)
#biological_analysis = biological_analysis[order(biological_analysis$Adjusted.P.value), ]
#biological_analysis = biological_analysis[1:15, ]
withr::with_dir(file.path(TARGET_dir, 'rib_genes'), {
  term = 'translation (GO:0006412)'
  target_genes = biological_analysis[biological_analysis$Term == term, 'Genes']
  term = stringr::str_replace(term, "/", "-")
  target_genes = stringr::str_split(target_genes, ";")[[1]]
  #target_genes = stringr::str_to_title(target_genes)
  norm_exp = monocle3::normalized_counts(cds)
  norm_exp = as.matrix(norm_exp)
  index_list = vector()
  for(gene in target_genes) { 
    index_list = c(index_list, which(tolower(gene) == tolower(rownames(norm_exp))))
  }
  
  norm_exp = norm_exp[index_list, ]
  intersecting_genes = intersect(rownames(norm_exp), ChIP_genes$symbol)
  intersecting_genes = c('rib', intersecting_genes)
  norm_exp = monocle3::normalized_counts(cds)
  norm_exp = as.matrix(norm_exp)
  norm_exp = norm_exp[intersecting_genes, ]
  # this will change 
  #norm_exp = norm_exp[apply(norm_exp, MARGIN = 1, FUN = max) > 1, ]
  pt = monocle3::pseudotime(cds)
  pt = data.frame(pseudotime = pt)
  plot_df = cbind(pt, t(norm_exp[, rownames(pt)]))
  smoothed_df = data.frame()
  for(gene in colnames(plot_df)) {
    if(gene == 'pseudotime') {
      next
    }
    else {
      yy = ksmooth(plot_df[, 'pseudotime'], plot_df[, gene], kernel="normal", bandwidth = 1.5, x.points=plot_df[, 'pseudotime'])
      if(nrow(smoothed_df) == 0) {
        smoothed_df = data.frame('pseudotime' = yy$x)
      }
      smoothed_df[, gene] = yy$y
    }
  }
  smoothed_df$pseudotime = NULL
  smoothed_df = t(smoothed_df)
  scaled_exp = t(scale(t(smoothed_df)))
  sorted_genes = names(sort(apply(scaled_exp, MARGIN = 1, FUN = which.max)))
  scaled_exp = scaled_exp[sorted_genes, ]
  pdf(file = paste0(term, "_dynamic_gene_heatmap.pdf"), height = 10)
  pheatmap(scaled_exp[sorted_genes, ], cluster_cols = FALSE, cluster_rows = FALSE)
  dev.off()
  
  convert_line_plot <- function(scaled_exp) {
    plot_df = data.frame()
    for(gene in rownames(scaled_exp)) {
      temp_plot = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                             scaled_exp = scaled_exp[gene, ], 
                             gene = gene)
      plot_df = rbind(plot_df, temp_plot)
    }
    
    plot_df$pseudotime = plot_df$pseudotime / max(plot_df$pseudotime)
    return(plot_df)
  }
  
  plot_df = convert_line_plot(scaled_exp)
  plot_df$Type = NA
  plot_df[plot_df$gene == 'rib', 'Type'] = 'rib'
  plot_df[plot_df$gene != 'rib', 'Type'] = 'targets of rib'
  p<-ggplot(plot_df, aes(x=pseudotime, y=scaled_exp, group=gene)) +
    xlab("pseudotime") + 
    ylab("smoothed and scaled expression") +
    ggtitle(paste0('genes in ', term)) +
    geom_line(aes(color=Type)) + theme_bw() 
  
  ggsave(paste0(term, "_dynamic_gene_line.png"), plot = p, width = 10, height = 8)
  
  # get the average 
  no_rib_scaled = scaled_exp[rownames(scaled_exp) != 'rib', ]
  plot_df = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                       scaled_exp = scaled_exp['rib', ], 
                       gene = 'rib')
  
  temp_df = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                       scaled_exp =   apply(no_rib_scaled, MARGIN = 2, FUN = mean), 
                       gene = 'targets of rib (average)')
  plot_df = rbind(plot_df, temp_df)
  p<-ggplot(plot_df, aes(x=pseudotime, y=scaled_exp, group=gene)) +
    xlab("pseudotime") + 
    ylab("smoothed and scaled expression") +
    ggtitle(paste0('genes in ', term)) +
    geom_line(aes(color=gene)) + theme_bw() 
  ggsave(paste0(term, "_dynamic_gene_line_avg.png"), plot = p, width = 10, height = 8)
  
})

##########################################
# CrebA genes dynamics with golgi-vescile 
dir.create(file.path(TARGET_dir, 'CrebA_genes'))

biological_analysis = read.csv(file.path(TARGET_dir, 'sig_GO_biological_late.csv'), row.names = 1)
#biological_analysis = biological_analysis[order(biological_analysis$Adjusted.P.value), ]
#biological_analysis = biological_analysis[1:15, ]
withr::with_dir(file.path(TARGET_dir, 'CrebA_genes'), {
  term = 'Golgi vesicle transport (GO:0048193)'
  target_genes = biological_analysis[biological_analysis$Term == term, 'Genes']
  term = stringr::str_replace(term, "/", "-")
  target_genes = stringr::str_split(target_genes, ";")[[1]]
  #target_genes = stringr::str_to_title(target_genes)
  norm_exp = monocle3::normalized_counts(cds)
  norm_exp = as.matrix(norm_exp)
  index_list = vector()
  for(gene in target_genes) { 
    index_list = c(index_list, which(tolower(gene) == tolower(rownames(norm_exp))))
  }
  
  norm_exp = norm_exp[index_list, ]
  intersecting_genes = rownames(norm_exp)
  intersecting_genes = c('CrebA', intersecting_genes)
  norm_exp = monocle3::normalized_counts(cds)
  norm_exp = as.matrix(norm_exp)
  norm_exp = norm_exp[intersecting_genes, ]
  # this will change 
  #norm_exp = norm_exp[apply(norm_exp, MARGIN = 1, FUN = max) > 1, ]
  pt = monocle3::pseudotime(cds)
  pt = data.frame(pseudotime = pt)
  plot_df = cbind(pt, t(norm_exp[, rownames(pt)]))
  smoothed_df = data.frame()
  for(gene in colnames(plot_df)) {
    if(gene == 'pseudotime') {
      next
    }
    else {
      yy = ksmooth(plot_df[, 'pseudotime'], plot_df[, gene], kernel="normal", bandwidth = 1.5, x.points=plot_df[, 'pseudotime'])
      if(nrow(smoothed_df) == 0) {
        smoothed_df = data.frame('pseudotime' = yy$x)
      }
      smoothed_df[, gene] = yy$y
    }
  }
  smoothed_df$pseudotime = NULL
  smoothed_df = t(smoothed_df)
  scaled_exp = t(scale(t(smoothed_df)))
  sorted_genes = names(sort(apply(scaled_exp, MARGIN = 1, FUN = which.max)))
  scaled_exp = scaled_exp[sorted_genes, ]
  pdf(file = paste0(term, "_dynamic_gene_heatmap.pdf"), height = 10)
  pheatmap(scaled_exp[sorted_genes, ], cluster_cols = FALSE, cluster_rows = FALSE)
  dev.off()
  
  convert_line_plot <- function(scaled_exp) {
    plot_df = data.frame()
    for(gene in rownames(scaled_exp)) {
      temp_plot = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                             scaled_exp = scaled_exp[gene, ], 
                             gene = gene)
      plot_df = rbind(plot_df, temp_plot)
    }
    
    plot_df$pseudotime = plot_df$pseudotime / max(plot_df$pseudotime)
    return(plot_df)
  }
  
  plot_df = convert_line_plot(scaled_exp)
  plot_df$Type = NA
  plot_df[plot_df$gene == 'CrebA', 'Type'] = 'CrebA'
  plot_df[plot_df$gene != 'CrebA', 'Type'] = 'vesicle transport related genes'
  p<-ggplot(plot_df, aes(x=pseudotime, y=scaled_exp, group=gene)) +
    xlab("pseudotime") + 
    ylab("smoothed and scaled expression") +
    ggtitle(paste0('genes in ', term)) +
    geom_line(aes(color=Type)) + theme_bw() 
  
  ggsave(paste0(term, "_dynamic_gene_line.png"), plot = p, width = 10, height = 8)
  
  # get the average 
  no_rib_scaled = scaled_exp[rownames(scaled_exp) != 'CrebA', ]
  plot_df = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                       scaled_exp = scaled_exp['CrebA', ], 
                       gene = 'CrebA')
  
  temp_df = data.frame(pseudotime = seq(1, ncol(scaled_exp)), 
                       scaled_exp =   apply(no_rib_scaled, MARGIN = 2, FUN = mean), 
                       gene = 'vesicle transport related genes (average)')
  plot_df = rbind(plot_df, temp_df)
  p<-ggplot(plot_df, aes(x=pseudotime, y=scaled_exp, group=gene)) +
    xlab("pseudotime") + 
    ylab("smoothed and scaled expression") +
    ggtitle(paste0('genes in ', term)) +
    geom_line(aes(color=gene)) + theme_bw() 
  ggsave(paste0(term, "_dynamic_gene_line_avg.png"), plot = p, width = 10, height = 8)
})
