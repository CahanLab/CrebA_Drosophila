import pandas as pd
import numpy as np 
import datetime
import os

# load in the scDE and maDE
sc_DE = pd.read_csv("../../analysis/results/v19/DE_genes_early_crebA_wt/Salivary Gland/mut_DE_genes.csv", index_col=0)
ma_DE = pd.read_excel("../input/CrebA_microarry_data/CrebA_all microarray data.xlsx")
promoter_tab = pd.read_csv("../output/tss_table/regulatory_regions.bed", sep = '\t', header = None)

# create the directory if needed 
output_path = "../output/find_bound_DE_genes/"
if os.path.isdir(output_path) == False:
    os.makedirs(output_path)

# get the DE genes 
reg_genes = dict()
for tmp_type in ['up', 'down']:
    if tmp_type == 'up':
        sub_ma_DE = ma_DE.loc[np.logical_and(ma_DE['CrebA vs. WT p-value'] < 0.05, ma_DE['CrebA vs. WT linear FC'] > 1.5), :]
        sub_sc_DE = sc_DE.loc[np.logical_and(sc_DE['logFC'] > 0.15, sc_DE['padj'] < 0.05), :]
    else:
        sub_ma_DE = ma_DE.loc[np.logical_and(ma_DE['CrebA vs. WT p-value'] < 0.05, ma_DE['CrebA vs. WT linear FC'] < -1.5), :]
        sub_sc_DE = sc_DE.loc[np.logical_and(sc_DE['logFC'] < -0.15, sc_DE['pval'] < 0.05), :]      
    sc_genes = np.array(sub_sc_DE['feature'])
    ma_genes = []

    # check if MA DE genes is in the transcriptome or not 
    sub_ma_DE['in_new_transcriptome'] = np.array(sub_ma_DE['Gene Symbol'].isin(np.array(promoter_tab[6])))
    sub_ma_DE.to_csv(os.path.join(output_path, tmp_type + '_ma_DE.csv'))
    for tmp_gene in list(sub_ma_DE['Gene Symbol']):
        if tmp_gene == datetime.datetime(2009, 9, 5, 0, 0): 
            ma_genes = ma_genes + ['Sep5']
        elif tmp_gene == datetime.datetime(2009, 9, 1, 0, 0): 
            ma_genes = ma_genes + ['Sep1']
        elif tmp_gene == datetime.datetime(2009, 9, 2, 0, 0): 
            ma_genes = ma_genes + ['Sep2']
        elif tmp_gene == datetime.datetime(2009, 12, 1, 0, 0): 
            ma_genes = ma_genes + ['dec-1']
        else:
            ma_genes = ma_genes + [x.strip() for x in tmp_gene.split("///")]
    ma_genes = np.array(ma_genes)
    ma_genes = ma_genes[ma_genes != '---']
    my_df = pd.DataFrame()
    my_df['genes'] = np.unique(np.concatenate((sc_genes, ma_genes), axis=None))
    my_df['MA_DE'] = my_df['genes'].isin(ma_genes)
    my_df['SC_DE'] = my_df['genes'].isin(sc_genes)
    my_df['both_DE'] = np.logical_and(my_df['MA_DE'], my_df['SC_DE'])
    my_df['in_transcriptome'] = np.array(my_df['genes'].isin(np.array(promoter_tab[6])))
    reg_genes[tmp_type] = my_df

# see if these genes are bound or not 
crebA_bound_df = pd.read_csv("../output/match_nearest_gene/fkh_sage_intersect_genes.csv", index_col = 0)
crebA_bound_df = crebA_bound_df.loc[crebA_bound_df['in_region_gene'].isna() == False, :]
bound_genes = list()
for i in crebA_bound_df.index: 
    bound_genes = bound_genes + crebA_bound_df.loc[i, 'in_region_gene'].split(",")

for tmp_type in ['up', 'down']:
    reg_genes[tmp_type]['bound'] = False
    reg_genes[tmp_type].loc[reg_genes[tmp_type]['genes'].isin(bound_genes), 'bound'] = True
    reg_genes[tmp_type].to_csv(os.path.join(output_path, tmp_type + "_DE.csv"))
