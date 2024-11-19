import pandas as pd
import numpy as np 
import datetime
import os

# load in the scDE and maDE
sc_DE = pd.read_csv("../../analysis/results/v19/DE_genes_early_crebA_wt/Salivary Gland/mut_DE_genes.csv", index_col=0)
ma_DE = pd.read_excel("../input/CrebA_microarry_data/CrebA_all microarray data.xlsx")
promoter_tab = pd.read_csv("../output/tss_table/promoter_regions.bed", sep = '\t', header = None)

# create the directory if needed 
output_path = "../output/find_bound_DE_genes/"
if os.path.isdir(output_path) == False:
    os.makedirs(output_path)

# get the DE genes 
# save the MA DE genes in the data first 
for tmp_type in ['up', 'down']:
    if tmp_type == 'up':
        sub_ma_DE = ma_DE.loc[np.logical_and(ma_DE['p-value(Genotype)'] < 0.05, ma_DE['CrebA vs. WT linear FC'] > 1.25), :] # might have to change the threshold a little bit 
        sub_sc_DE = sc_DE.loc[np.logical_and(sc_DE['logFC'] > 0.15, sc_DE['pval'] < 0.05), :]
    else:
        sub_ma_DE = ma_DE.loc[np.logical_and(ma_DE['p-value(Genotype)'] < 0.05, ma_DE['CrebA vs. WT linear FC'] < -1.25), :]
        sub_sc_DE = sc_DE.loc[np.logical_and(sc_DE['logFC'] < -0.15, sc_DE['pval'] < 0.05), :]      
    sc_genes = np.array(sub_sc_DE['feature'])
    ma_genes = []

    # check if MA DE genes is in the transcriptome or not 
    sub_ma_DE['in_new_transcriptome'] = np.array(sub_ma_DE['Gene Symbol'].isin(np.array(promoter_tab[6])))
    sub_ma_DE.to_csv(os.path.join(output_path, tmp_type + '_ma_DE.csv'))

