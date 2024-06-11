TARGET_dir = file.path("results", ANALYSIS_VERSION, "harmonized_wildtype_data")
dir.create(TARGET_dir)

##### reharmnoize the early data #####
wt_object = readRDS("accessory_data/wild_type_seurats/GSE234602_stage10_12_seurat_object.rds")
wt_object@meta.data$new_celltypes = wt_object@meta.data$manual_celltypes
wt_object@meta.data[wt_object@meta.data$manual_celltypes == 'Ant And Post Gut Primordium', 'new_celltypes'] = 'Ant and Post Gut Primordium'
wt_object@meta.data[wt_object@meta.data$manual_celltypes == 'Posterior Midgut Primordium', 'new_celltypes'] = 'Gut Endoderm'
wt_object@meta.data[wt_object@meta.data$manual_celltypes == 'Unknown (CNS)', 'new_celltypes'] = 'Unknown'
saveRDS(wt_object, file.path(TARGET_dir, 'stage10-12_reharmonized_seurat.rds'))

###### reharmonize the late data ######
wt_object = readRDS("accessory_data/wild_type_seurats/GSE234602_stage13_16_seurat_object.rds")
wt_object@meta.data$new_celltypes = wt_object@meta.data$manual_celltypes
wt_object@meta.data[wt_object@meta.data$manual_celltypes == 'Dorsal Epidermis', 'new_celltypes'] = 'Epidermis'
wt_object@meta.data[wt_object@meta.data$manual_celltypes == 'Ventral Epidermis', 'new_celltypes'] = 'Epidermis'
wt_object@meta.data[wt_object@meta.data$manual_celltypes == 'Late CNS', 'new_celltypes'] = 'CNS'
wt_object@meta.data[wt_object@meta.data$manual_celltypes == 'Early CNS', 'new_celltypes'] = 'CNS'
wt_object@meta.data[wt_object@meta.data$manual_celltypes == 'Lateral Sensory Neurons', 'new_celltypes'] = 'Sensory Neurons'
saveRDS(wt_object, file.path(TARGET_dir, 'stage13-16_reharmonized_seurat.rds'))
