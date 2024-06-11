library(Seurat)
library(harmony)
library(magrittr)
library(rvest)

SAMPLETYPE = "wt_rib"
TARGET_dir = file.path("results", ANALYSIS_VERSION, "wt_rib_integrated")
dir.create(TARGET_dir)

wt_rep1 = readRDS(file.path("results", ANALYSIS_VERSION,  "wt_rep1", "object.Rdata"))
wt_rep1@meta.data$batch = 'wt_rep1'

wt_rep2 = readRDS(file.path("results", ANALYSIS_VERSION,  "wt_rep2", "object.Rdata"))
wt_rep2@meta.data$batch = 'wt_rep2'

rib_rep1 = readRDS(file.path("results", ANALYSIS_VERSION,  "rib_rep1", "object.Rdata"))
rib_rep1@meta.data$batch = 'rib_rep1'

rib_rep2 = readRDS(file.path("results", ANALYSIS_VERSION,  "rib_rep2", "object.Rdata"))
rib_rep2@meta.data$batch = 'rib_rep2'

merged_seurat = merge(
  x = wt_rep1, 
  y = list(wt_rep2, rib_rep1, rib_rep2), 
  add.cell.ids = c("wt_rep1", 'wt_rep2', 'rib_rep1', 'rib_rep2')
)

raw_object = Seurat::CreateSeuratObject(merged_seurat@assays$RNA@counts, project = 'merged_samples')
raw_object@meta.data$batch = merged_seurat@meta.data$batch

saveRDS(raw_object, file = file.path(TARGET_dir, paste0(SAMPLETYPE, "_raw_merged.rds")))

object = Seurat::NormalizeData(raw_object)
cellCycleMarkers = read.csv("accessory_data/cellCycleMarkers.csv", skip = 1, header = T)
object %<>% CellCycleScoring(s.features = cellCycleMarkers$S.phase.markers., g2m.features = cellCycleMarkers$G2.M.phase.markers.)

genes_in_object = rownames(object@assays$RNA@meta.features)

mitochondrially_encoded_genes = read.table("accessory_data/mitochondrially_encoded_genes.tsv", header = F)[[1]]
ribosomal_rna_genes           = read.table("accessory_data/ribosomal_rna_genes.tsv", header = F)[[1]]
ribosomal_protein_genes       = grep("^rpl|^rps", ignore.case = T, value = T, genes_in_object)
object[["mitochondrial_transcript_total_expression"]] =
  GetAssayData(object, "counts") %>%
  extract( convert_fbgn_to_symbol( mitochondrially_encoded_genes ), ) %>%
  colSums
object[["ribosomal_rna_total_expression"]] =
  GetAssayData(object, "counts") %>%
  extract( convert_fbgn_to_symbol( ribosomal_rna_genes ), ) %>%
  colSums
object[["ribosomal_protein_total_expression"]] =
  GetAssayData(object, "counts") %>%
  extract( ribosomal_protein_genes, ) %>%
  colSums
object[["mitochondrial_transcript_norm_expression"]] = object[["mitochondrial_transcript_total_expression"]] / object$nCount_RNA
object[["ribosomal_rna_norm_expression"           ]] = object[["ribosomal_rna_total_expression"           ]] / object$nCount_RNA
object[["ribosomal_protein_norm_expression"       ]] = object[["ribosomal_protein_total_expression"       ]] / object$nCount_RNA
object[["log10_nCount_RNA"]] = object[["nCount_RNA"]] %>% log10

object %<>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
object %<>% ScaleData(features = VariableFeatures(object = object))
object %<>% RunPCA(features = VariableFeatures(object = object), npcs = 100)
pc_cutoff = RMTstat::qmp( 1, ndf=length(Cells(object)), pdim=length(VariableFeatures(object)), var=1)
singular_values = slot(Reductions(object, slot = "pca"), "stdev")
is_significant = singular_values^2 > pc_cutoff
num_pc = sum(is_significant)

# after getting the significant PCs, rerun PCA 
object %<>% RunPCA(features = VariableFeatures(object = object), npcs = num_pc)
object %<>% FindNeighbors(dims = 1:num_pc)
object %<>% RunUMAP(dim = 1:num_pc)
object %<>% FindClusters()
pdf(file = file.path(TARGET_dir, "naive_integration.pdf"))
DimPlot(object, group.by = "batch")
dev.off()

# harmony integration 
object %<>% harmony::RunHarmony("batch")
object %<>% FindNeighbors(reduction = "harmony", dims = 1:num_pc)
object %<>% RunUMAP(reduction = "harmony", dim = 1:num_pc)
object %<>% FindClusters(resolution = 0.95)
pdf(file = file.path(TARGET_dir, "harmony_integration.pdf"))
DimPlot(object, group.by = "batch")
dev.off()

saveRDS(object, file = file.path(TARGET_dir, 'object.rds'))

##################################
# look at the proportion study 
proportion_study <- function(object) { 
  plot_df = matrix(nrow = 4, ncol = 3)
  colnames(plot_df) = c('cluster_id', 'batch', 'proportion')
  plot_df = as.data.frame(plot_df)
  plot_df$cluster_id = 'full_data'
  plot_df$batch = names(table(object$batch))
  plot_df$proportion = table(object$batch) / length(object$batch)
  
  meta_tab = object@meta.data
  for(temp_cluster in unique(object$seurat_clusters)) { 
    temp_meta = meta_tab[meta_tab$seurat_clusters == temp_cluster, ]
    temp_plot_df = matrix(nrow = length(names(table(temp_meta$batch))), ncol = 3)
    colnames(temp_plot_df) = c('cluster_id', 'batch', 'proportion')
    temp_plot_df = as.data.frame(temp_plot_df)
    temp_plot_df$cluster_id = temp_cluster
    temp_plot_df$batch = names(table(temp_meta$batch))
    temp_plot_df$proportion = table(temp_meta$batch) / length(temp_meta$batch)
    
    plot_df = rbind(plot_df, temp_plot_df)
  }
  return(plot_df)
}
myProp = proportion_study(object)
p = ggplot(myProp) +
  geom_bar( aes(x=cluster_id, y=proportion, fill=batch), stat="identity") +
  theme_bw() + 
  coord_flip()

ggsave(file.path(TARGET_dir, "proportion_cluster.pdf"), plot = p, width = 7, height = 7.35)

Idents(object)=object[["seurat_clusters"]][[1]]
differentially_expressed = FindAllMarkers(object, test.use = "bimod")
withr::with_dir(
  file.path(TARGET_dir), 
  {
    
    write.csv(differentially_expressed, "marker_genes.csv")
    differentially_expressed %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n=20, wt = avg_log2FC) %>%
      dplyr::arrange(cluster) %>%
      write.csv("top_marker_genes_1-20_integrated.csv")
    
    MakeBasicPlots(object, reduction = "umap", plotname = "basic_stats_UMAP.pdf")
  }
)


####################################################
# web-scrape to get some tentative calling 
# load in the dictionary to convert the gene names to flybase gene name 
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

# here is a quick robot that will scrape the website for image 
BDGP_robot_image <- function(stage_name = 'stage13-16', gene_name = 'FBgn0003254') {
  url = paste0('https://insitu.fruitfly.org/cgi-bin/ex/report.pl?ftype=10&ftext=', gene_name)
  fly_data <- url %>%
    rvest::read_html() %>%
    rvest::html_nodes(xpath='/html/body/div[1]/div[2]/div[1]/table[2]') %>%
    rvest::html_nodes("tbody") %>% rvest::html_children()
  
  if(length(fly_data) == 0) { 
    return('None')
  }
  
  #link_string_button = paste0('<p><a class="btn btn-primary" data-toggle="collapse" href="#', gene_name, '" role="button" aria-expanded="false" aria-controls="', gene_name, '">Show Images</a></p>')
  link_string_button = '<summary>see images</summary>'
  link_string = ''
  
  for(temp_index in seq(1, length(fly_data))) { 
    cur_stage = fly_data[[temp_index]] %>% 
      html_nodes("th") %>% 
      html_text()
    if(cur_stage == stage_name) { 
      links = fly_data[[temp_index]] %>% 
        rvest::html_nodes("td") %>% 
        rvest::html_nodes('a') %>% 
        rvest::html_attr("href")
      links = links[grepl("cgi-bin", links) == FALSE]
      links = paste0('https://insitu.fruitfly.org', links)
      
      for(link in links) { 
        link_string = paste0(link_string, '<img src="', link, '" height=200></img>')
      }
    }
  }
  link_string = paste0('<p>', link_string, '</p>')
  link_string = paste0(link_string_button, link_string)
  link_string = paste0('<details>', link_string, "</details>")
  return(link_string)
}
# here is a quick robot that will scrape the website 
BDGP_robot <- function(stage_name = 'stage13-16', gene_name = 'FBgn0003254') { 
  url = paste0('https://insitu.fruitfly.org/cgi-bin/ex/report.pl?ftype=10&ftext=', gene_name)
  fly_data <- url %>%
    rvest::read_html() %>%
    rvest::html_nodes(xpath='/html/body/div[1]/div[2]/div[1]/table[2]') %>%
    rvest::html_nodes("tbody") %>% rvest::html_children()
  
  if(length(fly_data) == 0) { 
    return(vector())
  }
  
  for(temp_index in seq(1, length(fly_data))) { 
    cur_stage = fly_data[[temp_index]] %>% 
      html_nodes("th") %>% 
      html_text()
    if(cur_stage == stage_name) { 
      item = fly_data[[temp_index]] %>% 
        html_nodes("td") %>% 
        html_text()
      item_string = paste(item, collapse = '')
      item_string = stringr::str_remove_all(item_string, "\t")
      cell_types = stringr::str_split(item_string, "\n")[[1]]
      cell_types = cell_types[cell_types != ""]
      cell_types = trimws(cell_types)
    }
  }
  return(unique(cell_types))
}


withr::with_dir(
  file.path(TARGET_dir), 
  { 
    top_marker_genes = read.csv("top_marker_genes_1-20_integrated.csv")
    top_marker_genes$cell_type = NULL
    top_marker_genes$flybase_id = NULL
    for(temp_index in rownames(top_marker_genes)) { 
      gene_interest = top_marker_genes[temp_index, 'gene']
      flybase_interest = as.character(gene_converter[gene_interest])
      top_marker_genes[temp_index, 'flybase_id'] = flybase_interest
      BDGP_cellTypes = BDGP_robot(gene_name = flybase_interest)
      
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
      
      num_genes_list = sort(num_genes_list, decreasing = TRUE)
      score_list = sort(score_list, decreasing = TRUE)
      
      cluster_vector = c(cluster_vector, temp_cluster)
      if(num_genes_list[names(score_list)[1]] < 2) { 
        cell_type_vector = c(cell_type_vector, 'Unknown')
        num_genes_vector = c(num_genes_vector, 0)
        gene_scores_vector = c(gene_scores_vector, 0)
      } else{ 
        cell_type_vector = c(cell_type_vector, names(score_list)[1])
        num_genes_vector = c(num_genes_vector, num_genes_list[names(score_list)[1]])
        gene_scores_vector = c(gene_scores_vector, score_list[names(score_list)[1]])
      }
    }
    
    clusterAnnotation = data.frame(cluster = cluster_vector, 
                                   n_supporting_genes = num_genes_vector, 
                                   cellType_score = gene_scores_vector, 
                                   annotation = cell_type_vector)
    write.csv(clusterAnnotation, "tentativeCellTypes_BDGP.csv")
  }
)

# to add the image links for each 
withr::with_dir(
  file.path(TARGET_dir), 
  {
    top_marker_genes = read.csv('top_marker_genes_1-20_cellTypes.csv', row.names = 1)
    top_marker_genes$Images = NULL
    for(temp_index in rownames(top_marker_genes)) { 
      gene_interest = top_marker_genes[temp_index, 'gene']
      flybase_interest = as.character(gene_converter[gene_interest])
      top_marker_genes[temp_index, 'flybase_id'] = flybase_interest
      image_links = BDGP_robot_image(gene_name = flybase_interest)
      
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

object@meta.data$tentativeCellType = NULL
for(temp_cluster in unique(clusterAnnotation$cluster)) { 
  object@meta.data[object@meta.data$seurat_clusters == temp_cluster, 'tentativeCellType'] = clusterAnnotation[clusterAnnotation$cluster == temp_cluster, 'annotation']
}



DimPlot(object, group.by = "tentativeCellType", label = T, label.size = 5) +
  coord_fixed() + 
  ggtitle("Tentative cell type")
ggsave(file.path( TARGET_dir, "tentativeCellTypes_BDGP.png" ), width = 20, height = 7 )

saveRDS(object, file = file.path( TARGET_dir, 'BDGP_automated_annotation_object.rds'))

FeaturePlot(object, features = 'CG14756')
FeaturePlot(object, features = 'CG13159')
FeaturePlot(object, features = 'CG14453')

object@meta.data$smaller_tentativeCellType = NULL
for(i in rownames(metadata)) { 
  i = as.numeric(i)
  individual_object = readRDS(file.path("results", ANALYSIS_VERSION, metadata$sample[i], 'object.Rdata'))
  individual_object_meta = individual_object@meta.data
  rownames(individual_object_meta) = paste0(metadata$sample[i], "_", rownames(individual_object_meta))
  
  individual_clusterAnnotation = read.csv(file.path("results", ANALYSIS_VERSION, metadata$sample[i], 'tentativeCellTypes_BDGP.csv'), row.names = 1)
  individual_object_meta$tentativeCellType = NULL
  for(temp_cluster in unique(individual_clusterAnnotation$cluster)) { 
    individual_object_meta[individual_object_meta$seurat_clusters == temp_cluster, 'tentativeCellType'] = individual_clusterAnnotation[individual_clusterAnnotation$cluster == temp_cluster, 'annotation']
  }
  
  object@meta.data[rownames(individual_object_meta), 'smaller_tentativeCellType'] = individual_object_meta$tentativeCellType
}

DimPlot(object, group.by = 'smaller_tentativeCellType')
withr::with_dir(
  file.path(TARGET_dir), 
  { 
    
    curate_plot_df <- function(int_seurat, target_cell = 'all') { 
      UMAP_coord = int_seurat@reductions$umap@cell.embeddings
      UMAP_coord = as.data.frame(UMAP_coord)
      UMAP_coord$target_cells = NA
      UMAP_coord$batch = int_seurat$batch
      if(target_cell == 'all') { 
        UMAP_coord$target_cells = int_seurat$smaller_tentativeCellType
      } else{
        UMAP_coord[int_seurat$smaller_tentativeCellType == target_cell, 'target_cells'] = target_cell
        UMAP_coord[is.na(UMAP_coord$target_cells), 'target_cells'] = 'other'
      }
      return(UMAP_coord)
    }
    
    dir.create("overlayPlots_UMAPS")
    for(unique_celltype in unique(object@meta.data$smaller_tentativeCellType)) {
      plot_df = curate_plot_df(object, unique_celltype)
      plot_df[plot_df$target_cells != 'other', 'target_cells'] = paste0(plot_df[plot_df$target_cells != 'other', 'batch'], "_", plot_df[plot_df$target_cells != 'other', 'target_cells'])
      p = ggplot(plot_df) + 
        geom_point(data = subset(plot_df, target_cells == 'other'), aes(x=UMAP_1, y=UMAP_2, color = target_cells, alpha = 0.8))+
        geom_point(data = subset(plot_df, target_cells != 'other'), aes(x=UMAP_1, y=UMAP_2, color = target_cells, alpha = 0.8))+ 
        scale_color_manual(values=RColorBrewer::brewer.pal(12, 'Set3')[-c(2)]) + 
        ggtitle(unique_celltype)+
        theme_bw()+
        theme(legend.text = element_text(size = 12))
      
      ggsave(file.path('overlayPlots_UMAPS', paste0(unique_celltype, ".png")), plot = p, device = 'png', width = 12, height = 8)
    }  
  }
)

withr::with_dir(
  file.path(TARGET_dir), 
  {
    Seurat::VlnPlot(object, features = 'nCount_RNA', pt.size=0, y.max = 10000)
    Seurat::VlnPlot(object, features = 'nFeature_RNA', pt.size=0)
    Seurat::VlnPlot(object, features = 'log10_nCount_RNA', pt.size=0)
    Seurat::VlnPlot(object, features = 'G2M.Score', pt.size=0)
    Seurat::VlnPlot(object, features = 'ribosomal_rna_total_expression', pt.size=0)
    Seurat::VlnPlot(object, features = 'mitochondrial_transcript_total_expression', pt.size=0)
  }
)
