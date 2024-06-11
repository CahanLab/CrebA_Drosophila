library(rvest)
object=readRDS(file.path(RESULTS, "object.Rdata"))
Idents(object)=object[["seurat_clusters"]][[1]]

# Crank out cluster markers
withr::with_dir(
  file.path(RESULTS), 
  {
    if('top_marker_genes_1-40.csv' %in% list.files() == FALSE) { 
      differentially_expressed = FindAllMarkers(object, test.use = "bimod")
      write.csv(differentially_expressed, "marker_genes.csv")
      differentially_expressed %>%
        dplyr::group_by(cluster) %>%
        dplyr::top_n(n=20, wt = avg_log2FC) %>%
        dplyr::arrange(cluster) %>%
        write.csv("top_marker_genes_1-40.csv")      
    } else { 
      differentially_expressed = read.csv("top_marker_genes_1-40.csv", row.names = 1)
    }
  }
)

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

# NOTE these are legacy robot code that will scrape. To avoid pissing off the people at BDGP, 
#it's better to just download the database and rn
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


# add 
withr::with_dir(
  file.path(RESULTS), 
  {
    top_marker_genes = read.csv("top_marker_genes_1-40.csv")
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
      
      # remove the no ubiquitous 
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

# to add the image links for each 
withr::with_dir(
  file.path(RESULTS), 
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

object@meta.data$tentativeCellType = NULL
for(temp_cluster in unique(clusterAnnotation$cluster)) { 
  object@meta.data[object@meta.data$seurat_clusters == temp_cluster, 'tentativeCellType'] = clusterAnnotation[clusterAnnotation$cluster == temp_cluster, 'annotation']
}

# remove the embryonic status 
object@meta.data$tentativeCellType = stringr::str_remove(object@meta.data$tentativeCellType, "embryonic ")
object@meta.data$tentativeCellType = stringr::str_remove(object@meta.data$tentativeCellType, "embryonic/larval ")

DimPlot(object, group.by = "tentativeCellType", label = T, label.size = 5) +
  coord_fixed() + 
  ggtitle("Tentative cell type")
ggsave(file.path( RESULTS, "tentativeCellTypes_BDGP.png" ), width = 20, height = 7 )

saveRDS(object, file = file.path( RESULTS, 'BDGP_automated_annotation_object.rds'))
