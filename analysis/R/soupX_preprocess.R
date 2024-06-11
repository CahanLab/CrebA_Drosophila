library(SoupX)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "soupX_preprocess")
dir.create(TARGET_dir)

for(i in 1:5){ 
  dir.create(file.path(TARGET_dir, metadata[i, 'sample']))
  sc = SoupX::load10X(paste0("../quantification/", metadata[i, 'cellranger']))
  
  # Estimate rho
  sc = autoEstCont(sc, forceAccept = TRUE)
  
  # Clean the data
  out = adjustCounts(sc)
  
}





