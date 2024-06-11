wt_scRNA=readRDS(file.path("results", ANALYSIS_VERSION, "object_with_scores.Rdata"))

# Train a classifier based on the data from the Zinzen group. 
#
# Karaiskos, N., Wahle, P., Alles, J., Boltengagen, A., Ayoub, S., Kipar, C., ... & Zinzen, R. P. (2017). 
# The Drosophila embryo at single-cell transcriptome resolution. Science, 358(6360), 194-199.
#
karaiskos_et_al = read.table("accessory_data/karaiskos_et_al/GSE95025_high_quality_cells_digital_expression.txt.gz")
karaiskos_et_al %<>% extract(gene_identifiers$symbol,)
karaiskos_et_al[is.na(karaiskos_et_al)] = 0
karaiskos_et_al %<>% CreateSeuratObject()
karaiskos_et_al %<>% NormalizeData()
karaiskos_et_al[["karaiskos_et_al_labels"]] = karaiskos_et_al@assays$RNA@counts %>% colnames %>% gsub(".*_", "", .)
karaiskos_et_al[["barcode"]] = Cells(karaiskos_et_al)

# Use SingleCellNet for robust feature extraction
library(singleCellNet)
commonGenes = intersect(rownames(GetAssayData(karaiskos_et_al)), rownames(GetAssayData(wt_scRNA)))
length(commonGenes)
set.seed(0)
karaiskos_et_al_classifier = scn_train(stTrain = karaiskos_et_al[[c("karaiskos_et_al_labels", "barcode")]], 
                                       expTrain = GetAssayData(karaiskos_et_al)[commonGenes,], 
                                       nTopGenes = 10, 
                                       nRand = 2, 
                                       nTrees = 1000, 
                                       nTopGenePairs = 25, 
                                       dLevel = "karaiskos_et_al_labels", 
                                       colName_samp = "barcode")
classifier_results = scn_predict(cnProc=karaiskos_et_al_classifier[['cnProc']],
                                 expDat=GetAssayData(wt_scRNA)[commonGenes,],
                                 nrand = 2)
wt_scRNA[["karaiskos_et_al_labels"]] = 
  classifier_results[,Cells(wt_scRNA)] %>% 
  apply(2, function(x) names(which.max(x)))
withr::with_dir(file.path("results", ANALYSIS_VERSION), {
  MakeBasicPlots(object = wt_scRNA, 
                 reduction = "umap", 
                 plotname = "karaiskos_et_al_labels.pdf",
                 qc_features = "karaiskos_et_al_labels")
})
# Performance on the training set is good, so we are not *under*fitting.
karaiskos_et_al_classifier$cnProc$classifier$confusion
