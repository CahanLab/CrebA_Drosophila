TARGET_dir = file.path("results", ANALYSIS_VERSION, "ct_color_palettes")
dir.create(TARGET_dir, recursive = TRUE)

# this is to get all the CrebA data 
early_CrebA = readRDS(file.path('results', ANALYSIS_VERSION, 'manual_annotation_crebA_early/manual_celltype_object.rds'))
late_CrebA = readRDS(file.path('results', ANALYSIS_VERSION, 'manual_annotation_crebA/manual_celltype_object.rds'))

# get the wildtype data 
early_wt = readRDS(file.path('results', ANALYSIS_VERSION, 'harmonized_wildtype_data/stage10-12_reharmonized_seurat.rds'))
late_wt = readRDS(file.path('results', ANALYSIS_VERSION, 'harmonized_wildtype_data/stage13-16_reharmonized_seurat.rds'))

# get all the cell types 
all_cts_tot = unique(c(early_CrebA@meta.data$manual_celltypes, 
                   late_CrebA@meta.data$manual_celltypes, 
                   early_wt@meta.data$new_celltypes, 
                   late_wt@meta.data$new_celltypes))

all_cts = unique(c(early_CrebA@meta.data$manual_celltypes, 
                   late_CrebA@meta.data$manual_celltypes))

# I used glasbey Python package to generate the colour palette
# glasbey.extend_palette("Set2", palette_size=34, lightness_bounds=(30, 60))
color_palette = c('#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3', '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3', '#696100', '#8a2db2', '#006d6d', '#a23504', '#315dca', '#796971', '#1c9a31', '#b69200', '#b671f7', '#c2419e', '#828e69', '#189eae', '#ae714d', '#796daa', '#316d3d', '#0c92eb', '#209675', '#a68aaa', '#864571', '#795535', '#71868a', '#31c661', '#8282ff', '#86aa14', '#4d617d', '#db7500')
names(color_palette) = all_cts_tot

saveRDS(color_palette, file.path(TARGET_dir, 'ct_color_palette.rds'))
