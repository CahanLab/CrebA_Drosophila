library(Seurat)
library(harmony)
library(magrittr)
library(ggplot2)
library(stringr)

TARGET_dir = file.path("results", ANALYSIS_VERSION, "manual_annotation_early_wt12")
dir.create(TARGET_dir)

object = readRDS(file.path("results", ANALYSIS_VERSION, "early_wt12_integrated/BDGP_automated_annotation_object.rds"))

DimPlot(object)
DimPlot(object, group.by = 'Integrated_tentativeCellType', label = TRUE)


# classic SG genes 
VlnPlot(object, features = 'sage')
VlnPlot(object, features = 'pip')
FeaturePlot(object, features = 'toe')

# classic tracheal genes 
VlnPlot(object, features = 'trh')
FeaturePlot(object, features = 'trh')

# midline glia
VlnPlot(object, features = 'wrapper')
FeaturePlot(object, features = 'wrapper')

# plasmatocyte
VlnPlot(object, features = 'Fkbp14')
FeaturePlot(object, features = 'Fkbp14')

#glia
VlnPlot(object, features = 'repo')
FeaturePlot(object, features = 'repo')

# subset of CNS (voltage gated)
VlnPlot(object, features = 'nrv3')
FeaturePlot(object, features = 'nrv3')
VlnPlot(object, features = 'Sh')
FeaturePlot(object, features = 'Sh')

# crystal cells 
VlnPlot(object, features = 'PPO2')
FeaturePlot(object, features = 'PPO2')

# fat body 
VlnPlot(object, features = 'bgm')
FeaturePlot(object, features = 'bgm')
VlnPlot(object, features = 'CG6415')
FeaturePlot(object, features = 'CG6415')

#CNS 
VlnPlot(object, features = 'scrt', pt.size = 0)
FeaturePlot(object, features = 'scrt')
VlnPlot(object, features = 'Nrt', pt.size = 0)
FeaturePlot(object, features = 'Nrt')
VlnPlot(object, features = 'nerfin-1', pt.size = 0)
FeaturePlot(object, features = 'nerfin-1')
VlnPlot(object, features = 'elav', pt.size = 0)
FeaturePlot(object, features = 'elav')

# sensory 
VlnPlot(object, features = 'sv', pt.size = 0)
FeaturePlot(object, features = 'sv')
VlnPlot(object, features = 'sens', pt.size = 0)
FeaturePlot(object, features = 'sens')

# sensory - notch dependent 
VlnPlot(object, features = 'E(spl)m7-HLH', pt.size = 0)
FeaturePlot(object, features = 'E(spl)m7-HLH')
VlnPlot(object, features = 'Ocho', pt.size = 0)
FeaturePlot(object, features = 'Ocho')
VlnPlot(object, features = 'BobA', pt.size = 0)
FeaturePlot(object, features = 'BobA')

# epidermis 
VlnPlot(object, features = 'CG9628', pt.size = 0)
FeaturePlot(object, features = 'CG9628')

# more epidermis 
VlnPlot(object, features = 'Wdr62', pt.size = 0)
FeaturePlot(object, features = 'Wdr62')
VlnPlot(object, features = 'ft', pt.size = 0)
FeaturePlot(object, features = 'ft')

# more epidermis + foregut 
VlnPlot(object, features = 'Wdr62', pt.size = 0)
FeaturePlot(object, features = 'Wdr62')
VlnPlot(object, features = 'Fas3', pt.size = 0)
FeaturePlot(object, features = 'Fas3')

# hindgut 
VlnPlot(object, features = 'otp', pt.size = 0)
FeaturePlot(object, features = 'otp')
VlnPlot(object, features = 'Gmap', pt.size = 0)
FeaturePlot(object, features = 'Gmap')

# amnioserosa	
VlnPlot(object, features = 'CG12011', pt.size = 0)
FeaturePlot(object, features = 'CG12011')
VlnPlot(object, features = 'CG7997', pt.size = 0)
FeaturePlot(object, features = 'CG7997')

# Malpighian tubule
VlnPlot(object, features = 'Fst', pt.size = 0)
FeaturePlot(object, features = 'Fst')
VlnPlot(object, features = 'CG3699', pt.size = 0)
FeaturePlot(object, features = 'CG3699')

# gut 
VlnPlot(object, features = 'Ance', pt.size = 0)
FeaturePlot(object, features = 'Ance')
VlnPlot(object, features = 'DNaseII', pt.size = 0)
FeaturePlot(object, features = 'DNaseII')

# germ cells 
VlnPlot(object, features = 'BigH1', pt.size = 0)
FeaturePlot(object, features = 'BigH1')
VlnPlot(object, features = 'Sod1', pt.size = 0)
FeaturePlot(object, features = 'Sod1')
VlnPlot(object, features = 'bru1', pt.size = 0)
FeaturePlot(object, features = 'bru1')
VlnPlot(object, features = 'pgc', pt.size = 0)
FeaturePlot(object, features = 'pgc')

# yolk nuclei 
VlnPlot(object, features = 'CG5171', pt.size = 0)
FeaturePlot(object, features = 'CG5171')
VlnPlot(object, features = 'Oatp58Dc', pt.size = 0)
FeaturePlot(object, features = 'Oatp58Dc')

# heart muscle 
VlnPlot(object, features = 'Prm', pt.size = 0)
FeaturePlot(object, features = 'Prm')

# somatic muscle 
VlnPlot(object, features = 'Act57B', pt.size = 0)
FeaturePlot(object, features = 'Act57B')
VlnPlot(object, features = 'sug', pt.size = 0)
FeaturePlot(object, features = 'sug')

# longtitude visceral muslce 
VlnPlot(object, features = 'tey', pt.size = 0)
FeaturePlot(object, features = 'tey')
VlnPlot(object, features = 'beat-IIa', pt.size = 0)
FeaturePlot(object, features = 'beat-IIa')

manual_tab = read.csv(file.path(TARGET_dir, 'manualCellType.csv'))
object@meta.data$manual_celltypes = NULL

for(temp_cluster in unique(manual_tab$cluster)) { 
  object@meta.data[object@meta.data$seurat_clusters == temp_cluster, 'manual_celltypes'] = trimws(manual_tab[manual_tab$cluster == temp_cluster, 'annotation'])
}

p = DimPlot(object, group.by = 'manual_celltypes', label = TRUE)
ggsave(filename = file.path(TARGET_dir, "Dan's best guess.png"), plot = p, width = 10, height = 8)
saveRDS(object, file = file.path(TARGET_dir, "manual_celltype_object1.rds"))

