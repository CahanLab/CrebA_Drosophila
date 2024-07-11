TARGET_dir = file.path("results", ANALYSIS_VERSION, "Figures/CrebA_logFC_norm_exp")
dir.create(TARGET_dir)

logFC = read.csv(file.path("results/", ANALYSIS_VERSION, "/Figures/plot_spcg_logFC/stage10-12_logFC.csv"), row.names = 1)
norm_exp = read.csv(file.path("results/", ANALYSIS_VERSION, "/Figures/plot_spcg_wt/wt_early/spcg_scale_exp.csv"), row.names = 1)

logFC = logFC[logFC$feature == 'CrebA', ]
norm_exp = norm_exp[norm_exp$features.plot == 'CrebA', ]
intersect_ct = intersect(logFC$celltype, norm_exp$id)

rownames(logFC) = logFC$celltype
rownames(norm_exp) = norm_exp$id

logFC = logFC[intersect_ct, ]
norm_exp = norm_exp[intersect_ct, ]

plot_df = data.frame(row.names = intersect_ct, 
                     'cell_types' = intersect_ct, 
                     'logFC' = logFC$logFC, 
                     'norm_exp' = norm_exp$avg.exp)

p = ggplot(plot_df, aes(x = logFC, y = norm_exp)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(aes(label=ifelse(norm_exp>2,as.character(cell_types),'')),hjust=0,vjust=0) +
  xlab("CrebA logFC (mutant vs wildtype)") + 
  ylab("CrebA normalized exp (wildtype)") +
  ggtitle("CrebA normalized exp (wildtype) vs logFC (mutant vs wildtype)") +
  theme_cowplot()
ggsave(filename = file.path(TARGET_dir, 'logFC_exp_crebA.png'), plot = p, width = 9, height = 5)
