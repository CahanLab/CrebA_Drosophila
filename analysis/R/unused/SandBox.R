source("R/set_up_environment.R")
# This script is for interactive use for data exploration
setwd("/Users/DAndrewLab/Desktop/wt_scRNA/analysis/sandbox")
wt_scRNA=readRDS("../v2/initial_exploration.Rdata")
FeaturePlot(wt_scRNA, features = "rib", coord.fixed = T)
