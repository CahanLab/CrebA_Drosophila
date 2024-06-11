library(Seurat)
library(harmony)
library(magrittr)
library(ggplot2)
library(stringr)

SAMPLETYPE = "wt"
TARGET_dir = file.path("results", ANALYSIS_VERSION, "parse_wt_muscle")
dir.create(TARGET_dir)

object = readRDS(file.path("results", ANALYSIS_VERSION, 'manual_annotation_wt13/manual_celltype_object2.rds'))

sub_object = subset(object, cells = rownames(object@meta.data)[grepl("Somatic Muscle", object@meta.data$manual_celltypes)])

rm(list = 'object')
sub_object = Seurat::NormalizeData(sub_object)
cellCycleMarkers = read.csv("accessory_data/cellCycleMarkers.csv", skip = 1, header = T)
sub_object %<>% CellCycleScoring(s.features = cellCycleMarkers$S.phase.markers., g2m.features = cellCycleMarkers$G2.M.phase.markers.)

genes_in_object = rownames(sub_object@assays$RNA@meta.features)

mitochondrially_encoded_genes = read.table("accessory_data/mitochondrially_encoded_genes.tsv", header = F)[[1]]
ribosomal_rna_genes           = read.table("accessory_data/ribosomal_rna_genes.tsv", header = F)[[1]]
ribosomal_protein_genes       = grep("^rpl|^rps", ignore.case = T, value = T, genes_in_object)
sub_object[["mitochondrial_transcript_total_expression"]] =
  GetAssayData(sub_object, "counts") %>%
  extract( convert_fbgn_to_symbol( mitochondrially_encoded_genes ), ) %>%
  colSums
sub_object[["ribosomal_rna_total_expression"]] =
  GetAssayData(sub_object, "counts") %>%
  extract( convert_fbgn_to_symbol( ribosomal_rna_genes ), ) %>%
  colSums
sub_object[["ribosomal_protein_total_expression"]] =
  GetAssayData(sub_object, "counts") %>%
  extract( ribosomal_protein_genes, ) %>%
  colSums
sub_object[["mitochondrial_transcript_norm_expression"]] = sub_object[["mitochondrial_transcript_total_expression"]] / sub_object$nCount_RNA
sub_object[["ribosomal_rna_norm_expression"           ]] = sub_object[["ribosomal_rna_total_expression"           ]] / sub_object$nCount_RNA
sub_object[["ribosomal_protein_norm_expression"       ]] = sub_object[["ribosomal_protein_total_expression"       ]] / sub_object$nCount_RNA
sub_object[["log10_nCount_RNA"]] = sub_object[["nCount_RNA"]] %>% log10

sub_object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
sub_object %<>% ScaleData(features = VariableFeatures(object = sub_object))
sub_object %<>% RunPCA(features = VariableFeatures(object = sub_object), npcs = 100)
pc_cutoff = RMTstat::qmp( 1, ndf=length(Cells(sub_object)), pdim=length(VariableFeatures(sub_object)), var=1)
singular_values = slot(Reductions(sub_object, slot = "pca"), "stdev")
is_significant = singular_values^2 > pc_cutoff
num_pc = sum(is_significant)

# after getting the significant PCs, rerun PCA 
sub_object %<>% RunPCA(features = VariableFeatures(object = sub_object), npcs = num_pc)
sub_object %<>% FindNeighbors(dims = 1:num_pc)
sub_object %<>% RunUMAP(dim = 1:num_pc)
sub_object %<>% FindClusters()
pdf(file = file.path(TARGET_dir, "naive_integration.pdf"))
DimPlot(sub_object, group.by = "batch")
dev.off()

# harmony integration 
sub_object %<>% harmony::RunHarmony("batch")
sub_object %<>% FindNeighbors(reduction = "harmony", dims = 1:num_pc)
sub_object %<>% RunUMAP(reduction = "harmony", dim = 1:num_pc)
sub_object %<>% FindClusters(resolution = 0.6)
pdf(file = file.path(TARGET_dir, "harmony_integration.pdf"))
DimPlot(sub_object, group.by = "batch")
dev.off()

DimPlot(sub_object)

# Crank out cluster markers
withr::with_dir(
  file.path(TARGET_dir), 
  {
    if('top_marker_genes_1-20.csv' %in% list.files() == FALSE) { 
      differentially_expressed = FindAllMarkers(sub_object, test.use = "bimod")
      write.csv(differentially_expressed, "marker_genes.csv")
      differentially_expressed %>%
        dplyr::group_by(cluster) %>%
        dplyr::top_n(n=20, wt = avg_log2FC) %>%
        dplyr::arrange(cluster) %>%
        write.csv("top_marker_genes_1-20.csv")      
    } else { 
      differentially_expressed = read.csv("top_marker_genes_1-20.csv", row.names = 1)
    }
  }
)

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

BDGP_database_ct = read.csv("accessory_data/BDGP_marker_genes/insitu_annot.csv", header = FALSE)
colnames(BDGP_database_ct) = c("gene_name1", "gene_name2", 'fly_base_id', 'stage', 'cell_type')
BDGP_ct_assign <- function(BDGP_database_ct, stage_name = 'stage13-16', gene_name = 'FBgn0003254') { 
  stage_list = list()
  stage_list[['stage13-16']] = 6
  stage_list[['stage11-12']] = 5
  stage_list[['stage9-10']] = 4
  stage_list[['stage7-8']] = 3
  stage_list[['stage4-6']] = 2
  stage_list[['stage1-3']] = 1
  
  sub_BDGP = BDGP_database_ct[BDGP_database_ct$fly_base_id == gene_name, ]
  sub_BDGP = sub_BDGP[sub_BDGP$stage == stage_list[[stage_name]], ]
  
  if(nrow(sub_BDGP) > 0) {
    return(as.vector(sub_BDGP$cell_type))
  } else {
    return('None')
  }
}

BDGP_database_image = read.csv("accessory_data/BDGP_marker_genes/insitu_images.csv", header = FALSE)
BDGP_database_image = BDGP_database_image[, seq(1, 7)]
BDGP_image_assign <- function(BDGP_database_image, stage_name = 'stage13-16', gene_name = 'FBgn0003254') {
  
  colnames(BDGP_database_image) = c("gene_name1", "gene_name2", 'gene_name3',
                                    'fly_base_id', 'project_name', 'links', 
                                    'stage')
  stage_list = list()
  stage_list[['stage13-16']] = 6
  stage_list[['stage11-12']] = 5
  stage_list[['stage9-10']] = 4
  stage_list[['stage7-8']] = 3
  stage_list[['stage4-6']] = 2
  stage_list[['stage1-3']] = 1
  
  sub_BDGP_image = BDGP_database_image[BDGP_database_image$fly_base_id == gene_name, ]
  sub_BDGP_image = sub_BDGP_image[sub_BDGP_image$stage == stage_list[[stage_name]], ]
  
  link_string_button = '<summary>see images</summary>'
  link_string = ''
  
  if(nrow(sub_BDGP_image) > 0) { 
    links = paste0('https://insitu.fruitfly.org/insitu_image_storage/', sub_BDGP_image$links)
    for(link in links) { 
      link_string = paste0(link_string, '<img src="', link, '" height=200></img>')
    }
  }
  link_string = paste0('<p>', link_string, '</p>')
  link_string = paste0(link_string_button, link_string)
  link_string = paste0('<details>', link_string, "</details>")
  return(link_string)
}

withr::with_dir(
  file.path(TARGET_dir), 
  {
    top_marker_genes = read.csv("top_marker_genes_1-20.csv")
    top_marker_genes$cell_type = NULL
    top_marker_genes$flybase_id = NULL
    for(temp_index in rownames(top_marker_genes)) { 
      gene_interest = top_marker_genes[temp_index, 'gene']
      flybase_interest = as.character(gene_converter[gene_interest])
      top_marker_genes[temp_index, 'flybase_id'] = flybase_interest
      BDGP_cellTypes = BDGP_ct_assign(BDGP_database_ct, gene_name = flybase_interest)
      
      print(gene_interest)
      print(paste(BDGP_cellTypes, collapse = ";"))
      
      top_marker_genes[temp_index, 'cell_type'] = paste(BDGP_cellTypes, collapse = ";")
    }
    
    write.csv(top_marker_genes, "top_marker_genes_1-20_cellTypes.csv")
    
    cluster_vector = vector()
    cell_type_vector = vector()
    num_genes_vector = vector()
    gene_scores_vector = vector()
    
    for(temp_cluster in unique(top_marker_genes$cluster)) {
      temp_top_marker = top_marker_genes[top_marker_genes$cluster == temp_cluster, ]
      num_genes_list = list()
      score_list = list()
      
      for(temp_index in rownames(temp_top_marker)) { 
        log2Fold = temp_top_marker[temp_index, "avg_log2FC"]
        possible_cellTypes = stringr::str_split(temp_top_marker[temp_index, "cell_type"], ";")[[1]]
        if(possible_cellTypes[1] != "") { 
          for(possible_cellType in possible_cellTypes) { 
            if(possible_cellType %in% names(num_genes_list)) {
              num_genes_list[[possible_cellType]] = num_genes_list[[possible_cellType]] + 1
              score_list[[possible_cellType]] = score_list[[possible_cellType]] + log2Fold
            }
            else { 
              num_genes_list[[possible_cellType]] = 1
              score_list[[possible_cellType]] = log2Fold
            }
          }
        }
      }
      
      num_genes_list = unlist(num_genes_list)
      score_list = unlist(score_list)
      
      # remove the no staining 
      num_genes_list = num_genes_list[names(num_genes_list) != 'no staining']
      score_list = score_list[names(score_list) != 'no staining']
      
      # remove the no ubiquitous 
      num_genes_list = num_genes_list[names(num_genes_list) != 'ubiquitous']
      score_list = score_list[names(score_list) != 'ubiquitous']
      
      # remove the genes that are None 
      num_genes_list = num_genes_list[names(num_genes_list) != 'None']
      score_list = score_list[names(score_list) != 'None']
      
      num_genes_list = sort(num_genes_list, decreasing = TRUE)
      score_list = sort(score_list, decreasing = TRUE)
      
      cluster_vector = c(cluster_vector, temp_cluster)
      if(num_genes_list[names(score_list)[1]] < 2) { 
        cell_type_vector = c(cell_type_vector, 'Unknown')
        num_genes_vector = c(num_genes_vector, 0)
        gene_scores_vector = c(gene_scores_vector, 0)
      } else{ 
        if(score_list[names(score_list)[1]] < 5) {
          cell_type_vector = c(cell_type_vector, 'Unknown')
          num_genes_vector = c(num_genes_vector, num_genes_list[names(score_list)[1]])
          gene_scores_vector = c(gene_scores_vector, score_list[names(score_list)[1]])
        }
        else{
          cell_type_vector = c(cell_type_vector, names(score_list)[1])
          num_genes_vector = c(num_genes_vector, num_genes_list[names(score_list)[1]])
          gene_scores_vector = c(gene_scores_vector, score_list[names(score_list)[1]])
        }
      }
    }
    
    clusterAnnotation = data.frame(cluster = cluster_vector, 
                                   n_supporting_genes = num_genes_vector, 
                                   cellType_score = gene_scores_vector, 
                                   annotation = cell_type_vector)
    write.csv(clusterAnnotation, "tentativeCellTypes_BDGP.csv")
  }
)

withr::with_dir(
  file.path(TARGET_dir), 
  {
    top_marker_genes = read.csv('top_marker_genes_1-20_cellTypes.csv', row.names = 1)
    top_marker_genes$Images = NULL
    for(temp_index in rownames(top_marker_genes)) { 
      gene_interest = top_marker_genes[temp_index, 'gene']
      flybase_interest = as.character(gene_converter[gene_interest])
      top_marker_genes[temp_index, 'flybase_id'] = flybase_interest
      image_links = BDGP_image_assign(BDGP_database_image, gene_name = flybase_interest)
      
      print(gene_interest)
      
      top_marker_genes[temp_index, 'Images'] = image_links
    }
    top_marker_genes = subset(top_marker_genes, select = c('avg_log2FC', 'p_val_adj', 'cluster', 'gene', 'flybase_id', 'cell_type', 'Images'))
    top_marker_genes_df = DT::datatable(top_marker_genes, escape = FALSE, fillContainer = TRUE, filter = 'top', style = 'auto', width = '100px', height = '200px')
    doc <- htmltools::tagList(
      htmltools::div(top_marker_genes_df, style = "width=auto;height=auto")
    )
    doc[[1]]$children[[1]]$width = 'auto'
    doc[[1]]$children[[1]]$height = '800px'
    htmltools::save_html(html = doc, file = "widgets.html")
  }
)

sub_object@meta.data$finer_annotation = NULL
for(temp_cluster in unique(clusterAnnotation$cluster)) { 
  sub_object@meta.data[sub_object@meta.data$seurat_clusters == temp_cluster, 'finer_annotation'] = clusterAnnotation[clusterAnnotation$cluster == temp_cluster, 'annotation']
}

# remove the embryonic status 
sub_object@meta.data$finer_annotation = stringr::str_remove(sub_object@meta.data$finer_annotation, "embryonic ")
sub_object@meta.data$finer_annotation = stringr::str_remove(sub_object@meta.data$finer_annotation, "embryonic/larval ")

DimPlot(sub_object, group.by = "finer_annotation", label = T, label.size = 5) +
  coord_fixed() + 
  ggtitle("Finer cell type")
ggsave(file.path(TARGET_dir, "Finer_cellTypes.png" ), width = 20, height = 7 )

DimPlot(sub_object, group.by = "seurat_clusters", label = T, label.size = 5) +
  coord_fixed() + 
  ggtitle("seurat_clusters")
ggsave(file.path(TARGET_dir, "seurat_clusters.png" ), width = 8, height = 6)

saveRDS(sub_object, file = file.path(TARGET_dir, 'finer_annotation_object.rds'))

DimPlot(sub_object, group.by = 'Phase')
ggsave(file.path(TARGET_dir, "Phase.png" ), width = 8, height = 6 )

FeaturePlot(sub_object, features = 'log10_nCount_RNA')
ggsave(file.path(TARGET_dir, "log10_nCount_RNA.png" ), width = 20, height = 7 )

cytoTrace_meta = read.csv(file.path("results", ANALYSIS_VERSION, "wt13_cytoTrace/Somatic Muscle/", "cytotraced_sample_tab.csv"), row.names = 1)
cytoTrace_meta = cytoTrace_meta[rownames(sub_object@meta.data), ]

sub_object@meta.data$ct_pt = cytoTrace_meta$ct_pseudotime
FeaturePlot(sub_object, features = 'ct_pt')
ggsave(file.path(TARGET_dir, "ct_pt.png" ), width = 20, height = 7 )

