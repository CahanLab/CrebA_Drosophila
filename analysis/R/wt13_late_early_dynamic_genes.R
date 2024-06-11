library(tradeSeq)
library(amap)
library(dtwclust)
library(monocle3)
library(pheatmap)
library(SeuratWrappers)

subtype_folders = c("wt_late_early_salivary_gland", 
                    "wt_late_early_trachea", 
                    "wt_late_early_germ")

# calculate the statistics of dynamically expressed genes -- slow process
for(subtype_folder in subtype_folders) {
  print(subtype_folder)
  TARGET_dir = file.path("results", ANALYSIS_VERSION, subtype_folder)
  
  cds = readRDS(file.path(TARGET_dir, 'monocle3_no_batch_correct_object.rds'))
  expRaw = cds@assays@data$counts
  pt<-as.data.frame(pseudotime(cds))
  colnames(pt)<-"pseudotime"
  cw<-as.matrix(rep(1,nrow(pt)))
  rownames(cw)<-rownames(pt)
  ts<-tradeSeq::fitGAM(as.matrix(expRaw),pseudotime=as.matrix(pt),cellWeights=cw)
  
  saveRDS(ts, file = file.path(TARGET_dir, "tradeseq_fitgam_batch_corrected.rds"))
  
  ATres<-associationTest(ts)
  saveRDS(ATres, file = file.path(TARGET_dir, 'raw_associationTest_batch_corrected.rds'))

  ATres = ATres[!is.na(ATres$pvalue), ]
  ATres$adj_p = p.adjust(ATres$pvalue)
  ATres = ATres[ATres$adj_p < 0.05, ]
  saveRDS(ATres, file = file.path(TARGET_dir, "significant_associationTest_batch_corrected.rds"))
  
  startRes <- startVsEndTest(ts)
  startRes$adj_p = p.adjust(startRes$pvalue)
  saveRDS(startRes, file = file.path(TARGET_dir, 'raw_startvsendtest_batch_correct.rds'))  
}

# identify significantly dynamic genes based on the statistics calculated above 
for(subtype_folder in subtype_folders) {
  print(subtype_folder)
  TARGET_dir = file.path("results", ANALYSIS_VERSION, subtype_folder)
  
  ATres = readRDS(file.path(TARGET_dir, "significant_associationTest_batch_corrected.rds"))
  write.csv(ATres, file = file.path(TARGET_dir, 'significant_associationTest_batch_corrected.rds.csv'))
  
  startRes = readRDS(file.path(TARGET_dir, 'raw_startvsendtest_batch_correct.rds'))
  write.csv(startRes, file = file.path(TARGET_dir, 'raw_startvsendtest_batch_correct.csv'))
  
  startRes = startRes[startRes$adj_p < 0.05, ]
  write.csv(startRes, file = file.path(TARGET_dir, 'significant_startvsendtest_batch_correct.csv'))
}

# plot dynamically expressed genes and scale expression 
for(subtype_folder in subtype_folders) {
  print(subtype_folder)
  TARGET_dir = file.path("results", ANALYSIS_VERSION, subtype_folder)
  ATres = readRDS(file.path(TARGET_dir, "significant_associationTest_batch_corrected.rds"))
  ATres = ATres[ATres$meanLogFC >= 1, ]
  ATres$gene = rownames(ATres)
  
  cds = readRDS(file.path(TARGET_dir, 'monocle3_no_batch_correct_object.rds'))
  norm_exp = monocle3::normalized_counts(cds)
  norm_exp = as.matrix(norm_exp)
  norm_exp = norm_exp[rownames(ATres), ]
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
  pdf(file = file.path(TARGET_dir, "dynamic_gene_heatmap.pdf"))
  pheatmap(scaled_exp, cluster_cols = FALSE, cluster_rows = TRUE)
  dev.off()
  saveRDS(scaled_exp, file = file.path(TARGET_dir, "scaled_expression.rds"))
}


for(subtype_folder in subtype_folders) {
  print(subtype_folder)
  TARGET_dir = file.path("results", ANALYSIS_VERSION, subtype_folder)
  cds = readRDS(file.path(TARGET_dir, 'monocle3_no_batch_correct_object.rds'))
  
  cds = monocle3::cluster_cells(cds, resolution = 1e-2)
  plot_cells(cds, color_cells_by = "cluster", label_cell_groups = FALSE, cell_size = 1)
  marker_test_res <- top_markers(cds, group_cells_by="cluster", 
                                 reference_cells=1000, cores=8)
  
  seurat_object = exportCDS(cds, export_to = 'Seurat')
  as.Seurat(cds)
}

library(enrichR)
enrichR::setEnrichrSite("FlyEnrichr")

for(subtype_folder in subtype_folders) {
  print(subtype_folder)
  TARGET_dir = file.path("results", ANALYSIS_VERSION, subtype_folder)
  
  startRes = readRDS(file.path(TARGET_dir, 'raw_startvsendtest_batch_correct.rds'))
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
  
  biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
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
  
  biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
  write.csv(biological_analysis, file = file.path(TARGET_dir, 'sig_GO_biological_late.csv'))
}



TARGET_dir = file.path("results", ANALYSIS_VERSION, 'wt_late_early_trachea')
cds_tracheal = readRDS(file.path(TARGET_dir, 'monocle3_no_batch_correct_object.rds'))
ATres_tracheal = readRDS(file.path(TARGET_dir, 'significant_associationTest_batch_corrected.rds'))
plot_cells(cds_tracheal,
           genes=c("trh", 'rib', 'Osi6', 'Egfr', 'CG13272', 'Osi17', 'CG3376', 'N', 'CycB', 'CycA', 'Cdk1', 'Cdc27'),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE, cell_size = 2)

norm_exp = monocle3::normalized_counts(cds_tracheal)
norm_exp = as.matrix(norm_exp)
norm_exp = norm_exp[rownames(ATres_tracheal), ]
pt = monocle3::pseudotime(cds_tracheal)
pt = data.frame(pseudotime = pt)

my_correlation = Dist(norm_exp, method = "pearson", nbproc = 2, diag = FALSE, upper = FALSE)
clusters = hclust(my_correlation, method = 'average')
plot(clusters)
cluster_labels = cutree(clusters, k = 20)


TARGET_dir = file.path("results", ANALYSIS_VERSION, 'wt_late_early_salivary_gland')
cds_sg = readRDS(file.path(TARGET_dir, 'monocle3_no_batch_correct_object.rds'))
plot_cells(cds_sg,
           genes=c("eyg", "CrebA", 'rib', 'wbl', 'toe', 'sage', 'fkh', 'sano', 'pip', 'RpS17', 'RpL24', 'RpL41', 'RpL13A'),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE, cell_size = 2)



norm_exp = monocle3::normalized_counts(cds_sg)
norm_exp = as.matrix(norm_exp)
pt = monocle3::pseudotime(cds_sg)
pt = data.frame(pseudotime = pt)

#scaled_exp = t(scale(t(norm_exp)))
#scaled_exp = scaled_exp[, rownames(pt)[order(pt$pseudotime)]]
scaled_exp = norm_exp
plot_df = cbind(pt, t(scaled_exp[, rownames(pt)]))


# if they want to see the dynamics 
yy = ksmooth(plot_df[, 'pseudotime'], plot_df[, 'sage'], kernel="normal", bandwidth = 3, x.points=plot_df[, 'pseudotime'])
new_plot_df = data.frame("pt" = yy$x, 
                         'exp' = yy$y)
ggplot(new_plot_df, aes(x=pt, y=exp)) + geom_line()


#################################
library(tradeSeq)
cds = readRDS(file.path(TARGET_dir, 'monocle3_no_batch_correct_object.rds'))
expRaw = cds@assays@data$counts
pt<-as.data.frame(pseudotime(cds))
colnames(pt)<-"pseudotime"
cw<-as.matrix(rep(1,nrow(pt)))
rownames(cw)<-rownames(sampTab)
ts<-tradeSeq::fitGAM(as.matrix(expRaw),pseudotime=as.matrix(pt),cellWeights=cw)
ATres<-associationTest(ts)
saveRDS(ATres, file = file.path(TARGET_dir, 'raw_associationTest_no_batch_corrected.rds'))
ATres = readRDS(file.path(TARGET_dir, 'raw_associationTest_no_batch_corrected.rds'))

ATres = ATres[!is.na(ATres$pvalue), ]
ATres$adj_p = p.adjust(ATres$pvalue)
ATres = ATres[ATres$adj_p < 0.05, ]
ATres$genes = rownames(ATres)

norm_exp = seurat_object@assays$RNA@data
norm_exp = norm_exp[ATres$genes, rownames(pt)[order(pt$pseudotime)]]

scaled_exp = t(scale(t(norm_exp)))
scaled_exp = scaled_exp[, rownames(pt)[order(pt$pseudotime)]]

# smoothing the data 
grnKsmooth<-function(
    expDat,
    cells,
    BW = .25){
  
  BW = min(BW, (max(cells$pseudotime)-min(cells$pseudotime))/10)
  
  t1 = cells$pseudotime
  names(t1) = as.vector(cells$cell_name)
  sort(t1)
  expDat = expDat[,names(t1)]
  
  ans = matrix(0, nrow=nrow(expDat), ncol=ncol(expDat))
  for(i in seq(nrow(expDat))){
    yy = ksmooth(t1, expDat[i,], kernel="normal", bandwidth = BW, x.points=t1)
    # changed from uniform spacing (n.points=ncol(expDat)) to defining x coordinates
    ans[i,] = yy$y
  }
  rownames(ans) = rownames(expDat)
  colnames(ans) = colnames(expDat)
  ans
}

expSmoothed <- grnKsmooth(as.matrix(norm_exp), pt, BW=0.03)

ans = ans[, rownames(pt)[order(pt$pseudotime)]]
pheatmap::pheatmap(ans, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE)

plot_df = cbind(pt, t(scaled_exp[, rownames(pt)]))

yy = ksmooth(plot_df[, 'pseudotime'], plot_df[, 'CG3655'], kernel="normal", bandwidth = 0.25, x.points=plot_df[, 'pseudotime'])

new_plot_df = data.frame("pt" = yy$x, 
                         'exp' = yy$y)

ggplot(new_plot_df, aes(x=pt, y=exp)) + geom_point()
