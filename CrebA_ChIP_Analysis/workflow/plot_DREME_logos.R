library(dplyr)
library(ggseqlogo)
library(XML)
library(xml2)
library(cowplot)
library(ggplot2)

###### create the helper functions ######
get_motif <- function(target_motif) {
  target_motif_list = target_motif |> 
                      xml2::xml_find_all('pos') |> 
                      xml2::xml_attrs()
  prop_mat = matrix(data = NA, nrow = length(target_motif_list), ncol = 4)
  colnames(prop_mat) = c('A', 'C', 'G', 'T')
  for(i in seq(1, length(target_motif_list))) {
    prop_mat[i, ] = as.numeric(target_motif_list[[i]])
  }
  return(t(prop_mat))
}

curate_logos <- function(xml_object) { 
  motifs_obj = xml2::xml_find_all(xml_object, '//motif')
  output_list = list()
  for(i in seq(1, length(motifs_obj))) {
    target_motif = motifs_obj[[i]]
    DREME_id = xml_attr(target_motif, 'alt')
    eval = as.numeric(xml_attr(target_motif, 'evalue'))
    seq_id = xml_attr(target_motif, 'seq')
    output_list[[DREME_id]] = list()
    output_list[[DREME_id]][['eval']] = eval
    output_list[[DREME_id]][['seq_id']] = seq_id
    output_list[[DREME_id]][['prop_mat']] = get_motif(target_motif)
  }
  return(output_list)
}

##### make logo sequence down #####
args = commandArgs(trailingOnly = TRUE)
input_path = args[1]
output_path = args[2]

dreme_xml <- xml2::read_xml(file.path(input_path, "dreme.xml"))  # Replace with your XML file path
curated_logo_list = curate_logos(dreme_xml)
for(tmp_name in names(curated_logo_list)) { 
  p = ggseqlogo(curated_logo_list[[tmp_name]]$prop_mat) + 
    ggtitle(paste0(curated_logo_list[[tmp_name]]$seq_id, " motif eval: ", curated_logo_list[[tmp_name]]$eval)) + 
    cowplot::theme_cowplot()
  ggsave(filename = file.path(output_path, paste0(tmp_name, "_", curated_logo_list[[tmp_name]]$seq_id, ".png")), plot = p, width = 4, height = 2)
}
