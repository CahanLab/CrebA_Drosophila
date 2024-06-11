library(stringr)
library(enrichR)
library(openxlsx)
library(enrichR)

enrichR::setEnrichrSite('FlyEnrichr')

# load in the necessary code to cell type 
gtf_last_field = unique(read.table("../quantification/reference_genome_info/dmel-all-r6.33.gtf.gz", sep = "\t")[["V9"]])
gene_metadata = data.frame(
  flybase_id =
    gtf_last_field %>%
    strsplit(";") %>%
    sapply(extract2, 1) %>%
    gsub("^.*gene_id ", "", .),
  symbol =
    gtf_last_field %>%
    strsplit(";") %>%
    sapply(extract2, 2) %>%
    gsub("^.*gene_symbol ", "", .)
)
gene_metadata = gene_metadata %>% dplyr::distinct()
gene_converter = as.vector(gene_metadata$flybase_id)
names(gene_converter) = as.vector(gene_metadata$symbol)

TARGET_dir = file.path("results", ANALYSIS_VERSION, 'ChIP_overlap')
dir.create(TARGET_dir)

excel_path = file.path("accessory_data/ChIP_Targets/SG-and-Tracheal-Rib_IDR_20220928.xlsx")
page_names = openxlsx::getSheetNames(excel_path)

gene_info = read.csv2("accessory_data/FlyBase/best_gene_summary_fb_2022_05_Edits.tsv", sep = '\t')

ChIP_genes = openxlsx::read.xlsx(excel_path, sheet = page_names[1], colNames = TRUE)
DE_genes = read.csv(file.path("results", ANALYSIS_VERSION, "salivary_gland/DE_genes_wt_rib_sg_wilcox.csv"))
DE_genes = DE_genes[DE_genes$p_val_adj < 0.05, ]
DE_genes$FlyBase = gene_converter[DE_genes$X]

overlap_genes = DE_genes[DE_genes$FlyBase %in% ChIP_genes$FBgn, ]
overlap_genes$group = NA
overlap_genes[overlap_genes$avg_log2FC > 0, 'group'] = 'wt'
overlap_genes[overlap_genes$avg_log2FC < 0, 'group'] = 'rib'

temp_info = gene_info[gene_info$FBgn_ID %in% overlap_genes$FlyBase, ]
rownames(temp_info) = temp_info$FBgn_ID

overlap_genes$summary = NA
overlap_genes$summary = NA
for(gene in rownames(temp_info)) {
  overlap_genes[overlap_genes$FlyBase == gene, 'summary'] = temp_info[gene, 'Summary']
}

write.csv(overlap_genes, file = file.path(TARGET_dir, paste0("salivary_glands_overlaps_wilcox_genes.csv")))

# this is for the tracheal system 
ChIP_genes = openxlsx::read.xlsx(excel_path, sheet = page_names[2], colNames = TRUE)
DE_genes = read.csv(file.path("results", ANALYSIS_VERSION, "tracheal_system/DE_genes_wt_rib_ts_wilcox.csv"))
DE_genes = DE_genes[DE_genes$p_val_adj < 0.05, ]
DE_genes$FlyBase = gene_converter[DE_genes$X]

overlap_genes = DE_genes[DE_genes$FlyBase %in% ChIP_genes$FBgn, ]
overlap_genes$group = NA
overlap_genes[overlap_genes$avg_log2FC > 0, 'group'] = 'wt'
overlap_genes[overlap_genes$avg_log2FC < 0, 'group'] = 'rib'

temp_info = gene_info[gene_info$FBgn_ID %in% overlap_genes$FlyBase, ]
rownames(temp_info) = temp_info$FBgn_ID

overlap_genes$summary = NA
for(gene in rownames(temp_info)) {
  overlap_genes[overlap_genes$FlyBase == gene, 'summary'] = temp_info[gene, 'Summary']
}

write.csv(overlap_genes, file = file.path(TARGET_dir, paste0("tracheal_system_overlaps_wilcox_genes.csv")))
################################33
page_names = page_names[1:2]
find_overlaps <- function(DE_genes, ChIP_genes, ct_name = 'SG', type = 'wt') {
  overlap_genes = DE_genes[DE_genes$FlyBase %in% ChIP_genes$FBgn, ]
  outlap_genes = DE_genes[!DE_genes$FlyBase %in% ChIP_genes$FBgn, ]
  
  write.csv(overlap_genes, file = file.path(TARGET_dir, paste0(ct_name, "_", type, "_overlap_genes.csv")))
  write.csv(outlap_genes, file = file.path(TARGET_dir, paste0(ct_name, "_", type, "_outlap_genes.csv")))
  
  enrichment_results = enrichR::enrichr(
    genes = overlap_genes$X, 
    databases = c(
      "GO_Biological_Process_2018"
    )
  )
  biological_analysis = enrichment_results$GO_Biological_Process_2018
  biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
  write.csv(biological_analysis, file = file.path(TARGET_dir, paste0(ct_name, "_", type, "_overlap_GO.csv")))
  
  enrichment_results = enrichR::enrichr(
    genes = outlap_genes$X, 
    databases = c(
      "GO_Biological_Process_2018"
    )
  )
  biological_analysis = enrichment_results$GO_Biological_Process_2018
  biological_analysis = biological_analysis[biological_analysis$Adjusted.P.value < 0.05, ]
  write.csv(biological_analysis, file = file.path(TARGET_dir, paste0(ct_name, "_", type, "_outlap_GO.csv")))
}

# let's do SG first 
ChIP_genes = openxlsx::read.xlsx(excel_path, sheet = page_names[1], colNames = TRUE)

#SG wt
DE_genes = read.csv(file.path("results", ANALYSIS_VERSION, "salivary_gland/markers_wt_sg.csv"))
DE_genes$FlyBase = gene_converter[DE_genes$X]
find_overlaps(DE_genes, ChIP_genes, 'SG', 'wt')

# SG rib
DE_genes = read.csv(file.path("results", ANALYSIS_VERSION, "salivary_gland/markers_rib_sg.csv"))
DE_genes$FlyBase = gene_converter[DE_genes$X]
find_overlaps(DE_genes, ChIP_genes, 'SG', 'rib')

# TS wt
ChIP_genes = openxlsx::read.xlsx(excel_path, sheet = page_names[2], colNames = TRUE)

DE_genes = read.csv(file.path("results", ANALYSIS_VERSION, "tracheal_system/markers_wt_sg.csv"))
DE_genes$FlyBase = gene_converter[DE_genes$X]
find_overlaps(DE_genes, ChIP_genes, 'TS', 'wt')

# TS rib
DE_genes = read.csv(file.path("results", ANALYSIS_VERSION, "tracheal_system/markers_rib_sg.csv"))
DE_genes$FlyBase = gene_converter[DE_genes$X]
find_overlaps(DE_genes, ChIP_genes, 'TS', 'rib')
