library(Seurat)
library(stringr)
library(Matrix)
library(pheatmap)
library(enrichR)
library(RColorBrewer)
library(openxlsx)

cell_type = 'tracheal_system'
TARGET_dir = file.path("results", ANALYSIS_VERSION, cell_type)

object = readRDS(file.path(TARGET_dir, "ct_object.rds"))

dir.create(file.path(TARGET_dir, "time_genes"))
TARGET_dir = file.path(TARGET_dir, "time_genes")
excel_path = file.path("accessory_data/continuum_drosophila_embryonic_development/Andrew_Marker_Genes", paste0(cell_type, "_time_genes.xlsx"))
page_names = openxlsx::getSheetNames(excel_path)

marker_genes_list = list()
union_genes = c()
all_genes = c()
for(page_name in page_names) {
  print(page_name)
  time_name = stringr::str_split(page_name, " hr")[[1]][1]
  time_name = paste0(" hr", time_name)
  time_name = stringr::str_remove_all(time_name, " ")
  time_name = stringr::str_replace_all(time_name, "-", "_")
  marker_page = openxlsx::read.xlsx(excel_path, sheet = page_name, colNames = FALSE)
  colnames(marker_page) = c("gene", 
                            "p_val", 
                            "avg_log2FC", 
                            "pct.1", 
                            "pct.2", 
                            "p_val_adj", 
                            "cluster", 
                            "time")
  rownames(marker_page) = marker_page$gene
  marker_page = marker_page[intersect(marker_page$gene, rownames(object@assays$RNA@meta.features)), ]
  marker_genes_list[[time_name]] = marker_page
  
  if(length(union_genes) == 0) { 
    union_genes = marker_page$gene
  }
  else {
    union_genes = intersect(union_genes, marker_page$gene)
  }
  
  all_genes = c(all_genes, marker_page$gene)
  object <- AddModuleScore(
    object = object,
    features = list(marker_page$gene),
    name = time_name
  )
  
  colnames(object@meta.data) <-
    gsub(x = colnames(object@meta.data)
         , pattern = paste0(time_name,1)
         , replacement = time_name
    )
  
  p = FeaturePlot(object, features = time_name)
  ggsave(filename = file.path(TARGET_dir, paste0(time_name, "_UMAP.png")), plot = p, width = 8, height = 6)
  
  p = VlnPlot(object, features = time_name, group.by = 'batch', pt.size = 0)
  ggsave(filename = file.path(TARGET_dir, paste0(time_name, "_vln.png")), plot = p, width = 8, height = 6)
}

# plotout the timepoint specific genes 
unique_genes = names(table(all_genes)[table(all_genes) == 1])

for(time_frame in names(marker_genes_list)) {
  print(time_frame)
  
  marker_page = marker_genes_list[[time_frame]]
  marker_page = marker_page[marker_page$gene %in% unique_genes, ]
  
  if(nrow(marker_page) == 0) {
    next
  }
  else{
    object <- AddModuleScore(
      object = object,
      features = list(marker_page$gene),
      name = time_name
    )
    
    colnames(object@meta.data) <-
      gsub(x = colnames(object@meta.data)
           , pattern = paste0(time_name,1)
           , replacement = time_name
      )
  }
  
  p = FeaturePlot(object, features = time_name)
  ggsave(filename = file.path(TARGET_dir, paste0(time_name, "_specific_UMAP.png")), plot = p, width = 8, height = 6)
  
  p = VlnPlot(object, features = time_name, group.by = 'batch', pt.size = 0)
  ggsave(filename = file.path(TARGET_dir, paste0(time_name, "_specific_vln.png")), plot = p, width = 8, height = 6)
}

# load in the differentially expressed genes 
TARGET_dir = file.path("results", ANALYSIS_VERSION, cell_type)
DE_genes = read.csv(file.path(TARGET_dir, "DE_genes_wt_rib_ts.csv"))
DE_genes = DE_genes[DE_genes$p_val_adj < 0.05, ]

DE_genes = DE_genes[DE_genes$X %in% all_genes, ]
DE_genes = DE_genes[order(abs(DE_genes$avg_log2FC), decreasing = TRUE), ]

rib_markers = DE_genes[DE_genes$avg_log2FC < 0, 'X']
enrichment_results = enrichR::enrichr(
  genes = rib_markers, 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018", 
    'KEGG_2019'
  )
)
TARGET_dir = file.path(TARGET_dir, "time_genes")
biological_analysis = enrichment_results$GO_Biological_Process_2018
biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
write.csv(biological_analysis, file = file.path(TARGET_dir, 'rib_biological_GO.csv'))

wt_markers = DE_genes[DE_genes$avg_log2FC > 0, 'X']
enrichment_results = enrichR::enrichr(
  genes = wt_markers, 
  databases = c(
    "GO_Biological_Process_2018", 
    "GO_Molecular_Function_2018", 
    'KEGG_2019'
  )
)
biological_analysis = enrichment_results$GO_Biological_Process_2018
biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
write.csv(biological_analysis, file = file.path(TARGET_dir, 'wt_biological_GO.csv'))

# get the portion of ribosomal genes at each timepoint 
ribo_percent_list = list()
for(time_name in names(marker_genes_list)) {
  temp_genes = marker_genes_list[[time_name]]
  
  percentage = sum(grepl('Rp', unlist(temp_genes$gene))) / nrow(temp_genes)
  ribo_percent_list[[time_name]] = percentage
}

