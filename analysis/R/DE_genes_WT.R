# load in the spcg 
spcg_tab = readxl::read_excel("accessory_data/SPCG_files/SPCG List.xlsx")
spcg_tab = as.data.frame(spcg_tab)

# this is specifically running for V18 
wt_1_seurat = readRDS("results/v18/wt_rep1/object.Rdata")
wt_2_seurat = readRDS("results/v18/wt_rep2/object.Rdata")

wt_1_labels = read.csv("results/v18/Bianca_wt_cluster/wt1_clusters.csv")
colnames(wt_1_labels) = c('cluster_number', 'cluster_identification')

wt_2_labels = read.csv("results/v18/Bianca_wt_cluster/wt2_clusters.csv")
colnames(wt_2_labels) = c('cluster_number', 'cluster_identification')

# assigning the cell-types 
assign_ct <- function(seurat_obj, annotate_st) { 
  seurat_obj@meta.data$Bianca_CT = NA
  for(cluster in unique(annotate_st$cluster_number)) { 
    seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters == str_trim(cluster), 'Bianca_CT'] = annotate_st[annotate_st$cluster_number == cluster, 'cluster_identification']
  }
  return(seurat_obj)
}

wt_1_seurat = assign_ct(wt_1_seurat, wt_1_labels)
wt_1_seurat$batch = 'wt_1'

wt_2_seurat = assign_ct(wt_2_seurat, wt_2_labels)
wt_2_seurat$batch = 'wt_2'

# the command used to be 
# combined_seurat = merge(x = wt_1_seurat, y = wt_2_seurat, project = 'batch')
# but I think the data are not normalized properly that way. This way hopefully it's better 

combined_seurat_raw = merge(x = wt_1_seurat, y = wt_2_seurat, project = 'batch')
combined_seurat = Seurat::CreateSeuratObject(combined_seurat_raw@assays$RNA@counts)
combined_seurat@meta.data = combined_seurat_raw@meta.data
object = Seurat::NormalizeData(object)

common_ct_list = intersect(wt_1_seurat$Bianca_CT, wt_2_seurat$Bianca_CT)

diff_genes_list = list()
for(common_ct in common_ct_list) { 
  print(common_ct)  
  sub_object = subset(x = combined_seurat, cells = rownames(combined_seurat@meta.data)[combined_seurat@meta.data$Bianca_CT == common_ct])
  diff_genes = FindMarkers(sub_object, group.by = 'batch', ident.1 = 'wt_1', ident.2 = 'wt_2')
  diff_genes_list[[common_ct]] = diff_genes
}

TARGET_dir = file.path("results", ANALYSIS_VERSION, "DE_genes_wt")
dir.create(TARGET_dir)

# Functional enrichment analysis
library("enrichR")
enrichR::setEnrichrSite("FlyEnrichr")

withr::with_dir(
  file.path(TARGET_dir), 
  { 
    saveRDS(diff_genes_list, file = 'diff_genes_list.rds')
    for(ct in names(diff_genes_list)) { 
      DE_genes = diff_genes_list[[ct]]
      write.csv(DE_genes, file = paste0(ct, "_DE_genes.csv"))
    }
  }
)

# find enrichment for wt_1 
dir.create(file.path(TARGET_dir, 'wt_1'))
withr::with_dir(
  file.path(file.path(TARGET_dir, 'wt_1')), 
  { 
    for(ct in names(diff_genes_list)) { 
      DE_genes = diff_genes_list[[ct]]
      DE_genes = DE_genes[DE_genes$p_val_adj < 0.05, ]
      DE_genes = DE_genes[DE_genes$avg_log2FC > 0.5, ]
      enrichment_results = enrichR::enrichr(
        genes = rownames(DE_genes), 
        databases = c(
          "GO_Biological_Process_2018", 
          "GO_Molecular_Function_2018"
        )
      )
      biological_analysis = enrichment_results$GO_Biological_Process_2018
      molecular_analysis = enrichment_results$GO_Molecular_Function_2018
      
      write.csv(biological_analysis, file = paste0(ct, "_biological_process.csv"))
      write.csv(molecular_analysis, file = paste0(ct, "_molecular_function.csv"))
      
    }
  }
)

# find enrichment for wt_2
dir.create(file.path(TARGET_dir, 'wt_2'))
withr::with_dir(
  file.path(file.path(TARGET_dir, 'wt_2')), 
  { 
    for(ct in names(diff_genes_list)) { 
      DE_genes = diff_genes_list[[ct]]
      DE_genes = DE_genes[DE_genes$p_val_adj < 0.05, ]
      DE_genes = DE_genes[DE_genes$avg_log2FC < -0.5, ]
      enrichment_results = enrichR::enrichr(
        genes = rownames(DE_genes), 
        databases = c(
          "GO_Biological_Process_2018", 
          "GO_Molecular_Function_2018"
        )
      )
      biological_analysis = enrichment_results$GO_Biological_Process_2018
      molecular_analysis = enrichment_results$GO_Molecular_Function_2018
      
      write.csv(biological_analysis, file = paste0(ct, "_biological_process.csv"))
      write.csv(molecular_analysis, file = paste0(ct, "_molecular_function.csv"))
      
    }
  }
)
