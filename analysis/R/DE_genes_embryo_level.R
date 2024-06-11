TARGET_dir = file.path("results", ANALYSIS_VERSION, "DE_genes_embryo_level")
dir.create(TARGET_dir)

##### load in the spcg #####
spcg_tab = readxl::read_excel("accessory_data/SPCG_files/SPCG List.xlsx")
spcg_tab = as.data.frame(spcg_tab)

##### load in the flybase id converter ######
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

##### this is for the wild type (early) ######
# load in the wildtype for the early data 
wt_object = readRDS("accessory_data/wild_type_seurats/stage_10-12/manual_celltype_object1.rds")
mut_object = readRDS(file.path("results", ANALYSIS_VERSION, "manual_annotation_crebA_early/manual_celltype_object.rds"))

wt_object@meta.data$experimental_condition = 'Wt'
mut_object@meta.data$experimental_condition = 'mut'

combined_obj = merge(wt_object, mut_object)
object = Seurat::CreateSeuratObject(combined_obj@assays$RNA@counts, project = 'find_diff')
object@meta.data$experimental_condition = combined_obj@meta.data$experimental_condition
object@meta.data$batch = paste0(combined_obj@meta.data$experimental_condition, "_", combined_obj@meta.data$batch)

object = Seurat::NormalizeData(object)
rank_sum_test = presto::wilcoxauc(object, group_by = "experimental_condition")
rank_sum_test = rank_sum_test[rank_sum_test$pct_in > 10 | rank_sum_test$pct_out > 10, ]

write.csv(rank_sum_test, file.path(TARGET_dir, "early_wt_CrebA_embroy_DE_genes.csv"))

##### this is for the wild type (late) #####
wt_object = readRDS("accessory_data/wild_type_seurats/stage_13-16/manual_celltype_object4.rds")
mut_object = readRDS(file.path("results", ANALYSIS_VERSION, "crebA_integrated/BDGP_automated_annotation_object.rds"))

wt_object@meta.data$experimental_condition = 'Wt'
mut_object@meta.data$experimental_condition = 'mut'

combined_obj = merge(wt_object, mut_object)
object = Seurat::CreateSeuratObject(combined_obj@assays$RNA@counts, project = 'find_diff')
object@meta.data$experimental_condition = combined_obj@meta.data$experimental_condition
object@meta.data$batch = paste0(combined_obj@meta.data$experimental_condition, "_", combined_obj@meta.data$batch)

object = Seurat::NormalizeData(object)
rank_sum_test = presto::wilcoxauc(object, group_by = "experimental_condition")
rank_sum_test = rank_sum_test[rank_sum_test$pct_in > 10 | rank_sum_test$pct_out > 10, ]

write.csv(rank_sum_test, file.path(TARGET_dir, "late_wt_CrebA_embroy_DE_genes.csv"))

###### compare with the microarry data ######
MA_data = readxl::read_excel("accessory_data/CrebA_microarry_data/CrebA_all microarray data.xlsx")
MA_data = MA_data[MA_data$`Gene Symbol` != '---', ]
MA_data$`Gene Symbol` = stringr::str_split_fixed(MA_data$`Gene Symbol`, " ///", n = 2) %>% `[`(, 1)
MA_data = MA_data[!duplicated(MA_data$`Gene Symbol`), ]
MA_data$logFC = log2(abs(MA_data$`CrebA vs. WT linear FC`)) * sign(MA_data$`CrebA vs. WT linear FC`)
MA_data = data.frame(MA_data)
rownames(MA_data) = MA_data$Gene.Symbol

##### look at early embryo data #####
rank_sum_test = read.csv(file.path(TARGET_dir, "early_wt_CrebA_embroy_DE_genes.csv"), row.names = 1)
rank_sum_test = rank_sum_test[rank_sum_test$group == 'mut', ]
rownames(rank_sum_test) = rank_sum_test$feature
i_genes = intersect(rownames(rank_sum_test), rownames(MA_data))

sub_MA_data = MA_data[i_genes, ]
sub_rank_test = rank_sum_test[i_genes, ]

plot_df = data.frame(MA_genes = sub_MA_data$logFC, 
                     sc_genes = sub_rank_test$logFC)
rownames(plot_df) = gene_converter[i_genes]
plot_df = plot_df[rownames(plot_df) %in% spcg_tab$`Drosophila FBgn`, ]

# Fit a linear model
model <- lm(sc_genes ~ MA_genes, data = plot_df)

# Extract the coefficients
intercept <- coef(model)[1]
slope <- coef(model)[2]

# Calculate R-squared
r_squared <- summary(model)$r.squared

# Create a label for the equation and R-squared
label <- sprintf("y = %.2fx + %.2f\nR² = %.2f", slope, intercept, r_squared)

p = ggplot(data = plot_df, aes(x = MA_genes, y = sc_genes)) +
  geom_point() +
  theme_cowplot() +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Early single-cell vs MA",
       x = "MA DE Genes LogFC",
       y = "SC De Genes logFC") + 
  annotate("text", x = max(plot_df$MA_genes) * 0.3, y = max(plot_df$sc_genes) * 0.3, label = label, hjust = 1)
ggsave(file.path(TARGET_dir, 'MA_early_SC_comparisons.png'), plot = p, width = 8, height = 6)

##### this is for the late embryo data ######
rank_sum_test = read.csv(file.path(TARGET_dir, "late_wt_CrebA_embroy_DE_genes.csv"), row.names = 1)
rank_sum_test = rank_sum_test[rank_sum_test$group == 'mut', ]
rownames(rank_sum_test) = rank_sum_test$feature
i_genes = intersect(rownames(rank_sum_test), rownames(MA_data))

sub_MA_data = MA_data[i_genes, ]
sub_rank_test = rank_sum_test[i_genes, ]

plot_df = data.frame(MA_genes = sub_MA_data$logFC, 
                     sc_genes = sub_rank_test$logFC)
rownames(plot_df) = gene_converter[i_genes]
plot_df = plot_df[rownames(plot_df) %in% spcg_tab$`Drosophila FBgn`, ]

# Fit a linear model
model <- lm(sc_genes ~ MA_genes, data = plot_df)

# Extract the coefficients
intercept <- coef(model)[1]
slope <- coef(model)[2]

# Calculate R-squared
r_squared <- summary(model)$r.squared

# Create a label for the equation and R-squared
label <- sprintf("y = %.2fx + %.2f\nR² = %.2f", slope, intercept, r_squared)

p = ggplot(data = plot_df, aes(x = MA_genes, y = sc_genes)) +
  geom_point() +
  theme_cowplot() +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Late single-cell vs MA",
       x = "MA DE Genes LogFC",
       y = "SC De Genes logFC") + 
  annotate("text", x = max(plot_df$MA_genes) * 0.3, y = max(plot_df$sc_genes) * 0.3, label = label, hjust = 1)
ggsave(file.path(TARGET_dir, 'MA_late_SC_comparisons.png'), plot = p, width = 8, height = 6)

