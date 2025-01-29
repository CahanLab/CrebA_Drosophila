# This global variable is accessed throughout to decide where to save files
RESULTS = file.path("results", ANALYSIS_VERSION, SAMPLE)
dir.create(RESULTS, recursive = T, showWarnings = F)

# Set up working environment
library("Seurat")
library("ggplot2")
library("dplyr")
library("magrittr")
library("Matrix")
library('stringr')
library("harmony")
library("fgsea")
library("readxl")
library("cowplot")
library("scales")
library("enrichR")
library("openxlsx")
enrichR::setEnrichrSite('FlyEnrichr')

# This logs package versions etc; good for reproducibility
sink("sessioninfo.txt")
print(sessionInfo())
sink()

# This helps convert between two typed of gene ID: fbgn and gene symbol.
gene_identifiers =
  read.table("../quantification/reference_genome_info/dmel-all-r6.33.gtf.gz", sep = "\t")[["V9"]] %>%
  strsplit(";") %>%
  sapply(magrittr::extract, 1:2) %>%
  t %>%
  gsub("^.* ", "", .) %>%
  as.data.frame %>%
  (function(X) X[!duplicated(X), ]) %>%
  set_colnames(c("fbgn", "symbol"))
s2fbgn = setNames(gene_identifiers$fbgn, gene_identifiers$symbol)
fbgn2s = setNames(gene_identifiers$symbol, gene_identifiers$fbgn)
convert_symbol_to_fbgn = function(symbols){
  s2fbgn[symbols]
}
convert_fbgn_to_symbol = function(symbols){
  fbgn2s[symbols]
}
cat("Testing name conversion.\n")
print(convert_symbol_to_fbgn(c("rib", "sage", "sage", "sens")))
print(convert_fbgn_to_symbol(convert_symbol_to_fbgn(c("rib", "sage", "sage", "sens"))))

# Global variable to help establish unified colorscales throughout the paper
COLORSCALES = list(
  phase = scales::hue_pal()(3) %>% setNames(c("G1", "G2M", "S"))
)
# I often plot basics like PC1, assigned cc phase, and total counts. 
MakeBasicPlots = function(object, 
                          reduction = "umap",
                          plotname = "umap_qc.pdf",
                          qc_features = c(
                            "PC1",
                            "seurat_clusters",
                            "integratedClusters",
                            "annotation",
                            "shortAnnotation",
                            "sample", 
                            "genotype", 
                            "Phase",
                            "DoubletFinder_call",
                            "barcode_contains",
                            "log10_nCount_RNA", 
                            "nCount_RNA", 
                            "nFeature_RNA",
                            "fraction_unspliced",
                            "total_transcripts_unspliced",
                            "total_transcripts_spliced",
                            "ribosomal_rna_norm_expression",
                            "ribosomal_protein_norm_expression",
                            "mitochondrial_transcript_norm_expression",
                            "DoubletFinder_score", 
                            "rib", 
                            "sage"
                          )){
  pdf(plotname, width = 10, height = 8)
  for( f in qc_features ){
    if(f=="Phase"){
      DimPlot(object, group.by = "Phase", label = T, cols = COLORSCALES$phase) %>% print
    }
    try({
      if(is.numeric(object[[f]][[1]])){
        FeaturePlot(object, features = f, coord.fixed = T, reduction = reduction) %>% print
      } else {
        DimPlot(    object, group.by = f, label = T, reduction = reduction) %>% print
      }
    }, silent = T)
  }
  dev.off()
}

#' Easily fit a random forest classifier on a Seurat object.
#'
#' To speed up computations, these functions use a predefined set of 200 random projections for feature extraction.
#' The dimension is fixed at 17875 based on our reference genome. 
set.seed(0)
random_projections = matrix( rnorm( nrow(gene_identifiers) * 200 ), nrow = 200)
BuildClassifier = function(object, labels = "seurat_clusters"){
  if(ncol(random_projections) != nrow(GetAssayData(object))){
    stop("Expected 17875 genes. Was a different reference used?\n")
  }
  randomForest::randomForest( x = as.matrix( t( random_projections %*% GetAssayData(object) ) ), 
                              y = as.factor(as.character(object[[labels]][[1]])))
}
ApplyClassifier = function(query_object, classifier){
  if(ncol(random_projections) != nrow(GetAssayData(object))){
    stop("Expected 17875 genes. Was a different reference used?\n")
  }
  predict(classifier, newdata = as.matrix( t( random_projections %*% GetAssayData(query_object) ) ) )
}

#' Fill in missing values.
#' 
fillNA = function(x, filler = 0){
  x[is.na(x)] = filler
  x
}


#' Merge cluster annotations into a Seurat object.
#' 
#' @param ... Passed to merge
#' @param clusterAnnotations Dataframe to be merged with metadata.
#' @param by.x @param by.y How to match up the object (x) and the new annotations (y). As in base R merge(). 
#' @param namesToAdd What new fields will be present in the output.  
#' @param newNames If you want to take a field named "foo" in clusterAnnotations
#' and put its contents into a field named "bar" in the object, set namesToAdd = "foo"
#' and set newNames = "bar". 
#' 
AnnotateClusters = function(seuratObject, clusterAnnotations,          
                            by.y       = colnames(clusterAnnotations)[ 1],
                            by.x = by.y,
                            namesToAdd = colnames(clusterAnnotations)[-1],
                            newNames = namesToAdd,
                            ...){
  seuratObject[["barcode"]] = Cells(seuratObject)
  for(j in seq_along(namesToAdd)){
    seuratObject[[newNames[[j]]]] =
      seuratObject[[c(by.x, "barcode")]] %>% 
      merge( clusterAnnotations[c(by.y, namesToAdd[[j]])], 
             by.x = by.x,
             by.y = by.y,
             all.x = T, 
             all.y = F ) %>%
      tibble::column_to_rownames("barcode") %>%
      extract(Cells(seuratObject), namesToAdd[[j]])
  }
  seuratObject
}

#' Create a 3d umap html widget.
#'
#' @param object, @param umapkey Seurat object with 3d umap already present in \code{object[[umapkey]]}
#' @param group.by A categorical variable stored in \code{object@metadata}.
#'  
save3DUMAPCategorical = function(object, group.by, umapkey = "umap3D"){
  X = cbind(Embeddings(object[[umapkey]]), object@meta.data) 
  colnames(X)[1:3] = paste0("umap3d_", 1:3)
  my_formula = as.formula(paste0("~", group.by))
  p = 
    plotly::plot_ly(data = X,
                    x = ~umap3d_1, y = ~umap3d_2, z = ~umap3d_3,
                    color = my_formula,  opacity = .5, 
                    colors = scales::hue_pal()(length(unique(X[[group.by]]))),
                    type = "scatter3d",
                    mode = "markers",
                    marker = list(size = 5, width=2),
                    hoverinfo="text") %>%
    plotly::add_trace(data = X %>%
                        dplyr::group_by_(group.by) %>%
                        dplyr::summarise(umap3d_1 = mean(umap3d_1),
                                         umap3d_2 = mean(umap3d_2),
                                         umap3d_3 = mean(umap3d_3)), 
                      x = ~umap3d_1, y = ~umap3d_2, z = ~umap3d_3, 
                      text = my_formula,
                      opacity = 1, color = ~1, colors = colorRamp("black"),
                      type = "scatter3d", 
                      mode = "markers+text",
                      marker = list(size = 10, width=4), 
                      textfont = list(color = '#000000', size = 16),
                      hoverinfo="text") 
  htmlwidgets::saveWidget(p, file = paste0(group.by, ".html"))
  return(p)
}

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

modified_dotPlot_df <- function(
    object,
    assay = NULL,
    features,
    cols = c("lightgrey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 6,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = 'radius',
    scale.min = NA,
    scale.max = NA
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% rownames(x = brewer.pal.info))
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(
      X = 1:length(features),
      FUN = function(x) {
        return(rep(x = names(x = features)[x], each = length(features[[x]])))
      }
    ))
    if (any(is.na(x = feature.groups))) {
      warning(
        "Some feature groups are unnamed.",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  
  data.features <- FetchData(object = object, vars = features, cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  } else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(
    X = unique(x = data.features$id),
    FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(
        X = data.use,
        MARGIN = 2,
        FUN = function(x) {
          return(mean(x = expm1(x = x)))
        }
      )
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    }
  )
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(
      what = rbind,
      args = lapply(X = data.plot, FUN = unlist)
    )
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(
    X = names(x = data.plot),
    FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    }
  )
  data.plot <- do.call(what = 'rbind', args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning(
      "Only one identity present, the expression values will be not scaled",
      call. = FALSE,
      immediate. = TRUE
    )
  } else if (ngroup < 5 & scale) {
    warning(
      "Scaling data with a low number of groups may produce misleading results",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  avg.exp.scaled <- sapply(
    X = unique(x = data.plot$features.plot),
    FUN = function(x) {
      data.use <- data.plot[data.plot$features.plot == x, 'avg.exp']
      if (scale) {
        data.use <- scale(x = data.use)
        data.use <- MinMax(data = data.use, min = col.min, max = col.max)
      } else {
        data.use <- log1p(x = data.use)
      }
      return(data.use)
    }
  )
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(
    x = data.plot$features.plot,
    levels = features
  )
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(
      X = as.character(x = data.plot$id),
      FUN = gsub,
      FUN.VALUE = character(length = 1L),
      pattern =  paste0(
        '^((',
        paste(sort(x = levels(x = object), decreasing = TRUE), collapse = '|'),
        ')_)'
      ),
      replacement = '',
      USE.NAMES = FALSE
    )
    data.plot$colors <- mapply(
      FUN = function(color, value) {
        return(colorRampPalette(colors = c('grey', color))(20)[value])
      },
      color = cols[splits.use],
      value = avg.exp.scaled
    )
  }
  color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, 'pct.exp'] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, 'pct.exp'] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(
      x = feature.groups[data.plot$features.plot],
      levels = unique(x = feature.groups)
    )
  }
  return(data.plot)
}

