library(singleCellNet)

seurat_obj = readRDS("accessory_data/annotations/2022JUNE23/Manual_wt_rep1_labelled/object.Rdata")
RESULTS = file.path("results", ANALYSIS_VERSION, 'SCN_training')
dir.create(RESULTS)

expression_mat = seurat_obj@assays$RNA@counts
expression_mat = as.matrix(expression_mat)

meta_tab = seurat_obj@meta.data
table(meta_tab$annotation)
table(meta_tab$shortAnnotation)

# this is for training 
withr::with_dir(
  file.path(RESULTS), 
  {
    set.seed(100)
    stList = splitCommon(sampTab = meta_tab, ncells = 100, dLevel = 'annotation')
    stTrain = stList[[1]]
    expTrain = expression_mat[, rownames(stTrain)]
    saveRDS(stList, file = 'stList.rds')
    
    system.time(class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 10, nRand = 70, nTrees = 1000, nTopGenePairs = 25, dLevel = "annotation", colName_samp = "barcode"))
    saveRDS(class_info, file = 'class_info.rds')
  }
)

withr::with_dir(
  file.path(RESULTS), 
  {
    #validate data
    stTestList = splitCommon(sampTab=stList[[2]], ncells=100, dLevel="annotation") #normalize validation data so that the assessment is as fair as possible
    stTest = stTestList[[1]]
    expTest = expression_mat[ ,rownames(stTest)]
    
    #predict
    classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 50)
    tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "barcode", classTrain = "annotation", classQuery = "annotation", nRand = 50)
    plot_PRs(tm_heldoutassessment)
    
  }
)
