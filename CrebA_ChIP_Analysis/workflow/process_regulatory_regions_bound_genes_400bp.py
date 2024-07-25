import pandas as pd
import numpy as np 
import os 

output_path = '../output/bound_genes_regulatory_regions_400bp'
os.makedirs(output_path, exist_ok=True)

# load in the data 
reg_region_tab = pd.read_csv("../output/tss_table/regulatory_regions_400bp.bed", sep = '\t', header=None)
all_bound_genes = pd.read_csv("../output/plot_location_crebA_binding/all_bound_genes_dist.csv", index_col=0)
target_bound_genes = all_bound_genes.loc[np.abs(all_bound_genes['dist_tss']) < 200, :].copy()

# load in the genes that are bound genes 
down_bound_genes_df = pd.read_csv("../output/find_bound_DE_genes/down_DE.csv", index_col = 0)
up_bound_genes_df = pd.read_csv("../output/find_bound_DE_genes/up_DE.csv", index_col = 0)

down_bound_genes_df = down_bound_genes_df.loc[down_bound_genes_df['SC_DE'] == True, :].copy()
down_bound_genes_df = down_bound_genes_df.loc[down_bound_genes_df['bound'] == True, :].copy()
up_bound_genes_df = up_bound_genes_df.loc[up_bound_genes_df['SC_DE'] == True, :].copy()
up_bound_genes_df = up_bound_genes_df.loc[up_bound_genes_df['bound'] == True, :].copy()

down_bound_genes_df = down_bound_genes_df.loc[down_bound_genes_df['genes'].isin(target_bound_genes['bound_gene']), :].copy()
up_bound_genes_df = up_bound_genes_df.loc[up_bound_genes_df['genes'].isin(target_bound_genes['bound_gene']), :].copy()

# make the down bed files 
sub_reg_tab = reg_region_tab.loc[reg_region_tab[6].isin(down_bound_genes_df['genes']), :]
reg_tab_bed = sub_reg_tab.loc[:, [0, 2, 3, 6]].copy()
reg_tab_bed.to_csv(os.path.join(output_path, 'down_genes_reg_regions.bed'), sep = '\t', header = None, index = None)

# make the down bed files 
sub_reg_tab = reg_region_tab.loc[reg_region_tab[6].isin(up_bound_genes_df['genes']), :]
reg_tab_bed = sub_reg_tab.loc[:, [0, 2, 3, 6]].copy()
reg_tab_bed.to_csv(os.path.join(output_path, 'up_genes_reg_regions.bed'), sep = '\t', header = None, index = None)

##### get the most static but also bound genes #####
sc_DE = pd.read_csv("../../analysis/results/v19/DE_genes_early_crebA_wt/Salivary Gland/mut_DE_genes.csv", index_col=0)
sc_DE = sc_DE.sort_values(by='pval', ascending = False)
                          
crebA_bound_df = pd.read_csv("../output/match_nearest_gene/fkh_sage_intersect_genes.csv", index_col = 0)
crebA_bound_df = crebA_bound_df.loc[crebA_bound_df['in_region_gene'].isna() == False, :]
bound_genes = list()
for i in crebA_bound_df.index: 
    bound_genes = bound_genes + crebA_bound_df.loc[i, 'in_region_gene'].split(",")

sc_DE['bound'] = False
sc_DE.loc[sc_DE['feature'].isin(bound_genes), 'bound'] = True
sc_DE = sc_DE.loc[sc_DE['bound'] == True, :]
sc_DE = sc_DE.loc[np.logical_or(sc_DE['pct_in'] > 30, sc_DE['pct_out'] > 30), :]

sc_DE = sc_DE.sort_values(by = 'pval', ascending=False)
sub_sc_DE = sc_DE.iloc[0:150, :].copy()

sub_reg_tab = reg_region_tab.loc[reg_region_tab[6].isin(sub_sc_DE['feature']), :]
reg_tab_bed = sub_reg_tab.loc[:, [0, 2, 3, 6]].copy()
reg_tab_bed.to_csv(os.path.join(output_path, 'static_genes_reg_regions.bed'), sep = '\t', header = None, index = None)

sc_DE.to_csv(os.path.join(output_path, 'static_DE.csv'))