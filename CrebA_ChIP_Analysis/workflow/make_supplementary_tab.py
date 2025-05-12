import pandas as pd
import os 
import numpy as np 

# this is to generate the supplementary table for the CrebA SG network 
DE_genes_path = "../output/find_bound_SG_DE_genes"
bound_genes_path = "../output/match_nearest_gene/fkh_sage_intersect_genes_1000.csv"


# get the spcgs 
spcg_df = pd.read_excel('../input/SPCG_files/SPCG List.xlsx')
conversion_tab = pd.read_csv("../input/flybase_gene_conversion/conversion_tab.csv")
conversion_tab.index = conversion_tab['flybase']
sub_conversion_tab = conversion_tab.loc[spcg_df['Drosophila FBgn'], :]
spcg_df['new_gene_id'] = np.array(sub_conversion_tab['gene_names'])
spcgs_genes = np.unique(spcg_df['new_gene_id'])

# get all the bound genes 
crebA_bound_df = pd.read_csv(bound_genes_path, index_col = 0)
concatenated_array = np.concatenate((np.array(crebA_bound_df['nearest_gene_1'].dropna()), np.array(crebA_bound_df['nearest_gene_2'].dropna())))
concatenated_array = np.unique(concatenated_array)
concatenated_array = list(concatenated_array)
for element in concatenated_array:
    # Split the string by whitespace
    split_elements = element.split(",")
    if len(split_elements) > 1:
        concatenated_array.extend(split_elements)
        concatenated_array.remove(element)
bound_genes = np.unique(concatenated_array)

# look at the all the other genes 
down_DE_genes_df = pd.read_csv(os.path.join(DE_genes_path, "down_DE.csv"), index_col=0)
down_genes = np.unique(down_DE_genes_df['genes'])

sc_down_DE_genes_df = down_DE_genes_df.loc[down_DE_genes_df['SC_DE'] == True, :]
sc_down_genes = np.unique(sc_down_DE_genes_df['genes'])

ma_down_DE_genes_df = down_DE_genes_df.loc[down_DE_genes_df['MA_DE'] == True, :]
MA_down_genes = np.unique(ma_down_DE_genes_df['genes'])

insitu_down_DE_genes_df = down_DE_genes_df.loc[down_DE_genes_df['in_situ_DE'] == True, :]
insitu_genes = np.unique(insitu_down_DE_genes_df['genes'])

# get the genes that are up 
up_DE_genes_df = pd.read_csv(os.path.join(DE_genes_path, "up_DE.csv"), index_col=0)
up_genes = np.unique(up_DE_genes_df['genes'])

sc_up_DE_genes_df = up_DE_genes_df.loc[up_DE_genes_df['SC_DE'] == True, :]
sc_up_genes = np.unique(sc_up_DE_genes_df['genes'])

ma_up_DE_genes_df = up_DE_genes_df.loc[up_DE_genes_df['MA_DE'] == True, :]
MA_up_genes = np.unique(ma_up_DE_genes_df['genes'])

# Concatenate multiple numpy arrays into one
all_genes = np.concatenate((sc_down_genes, MA_down_genes, insitu_genes, sc_up_genes, MA_up_genes, bound_genes))
all_genes = np.unique(all_genes)

SG_genes = pd.read_csv("../../analysis/results/v19/early_wt_gsea/Salivary Gland/markers_genes.csv", index_col = 0)
SG_genes = SG_genes.loc[SG_genes['pct.1'] >= 0.1, :]

mut_diff_genes = pd.read_csv("../../analysis/results/v19/DE_genes_early_crebA_wt/Salivary Gland/mut_DE_genes.csv", index_col = 0)
mut_diff_genes = mut_diff_genes.loc[mut_diff_genes['logFC'] < 0, :]

MA_SG = np.intersect1d(MA_down_genes, SG_genes.index)
MA_SG = np.intersect1d(MA_SG, mut_diff_genes['feature'])

# create a dataframe
supp_tab = pd.DataFrame(all_genes, columns = ['genes'])

# add the bound column
supp_tab['bound'] = supp_tab['genes'].isin(bound_genes)

# add the spcg column
supp_tab['SPCG'] = supp_tab['genes'].isin(spcgs_genes)

# add the activated genes 
activated_genes = np.unique(list(sc_down_genes) + list(insitu_genes) + list(MA_SG))

# add the repressed genes 
repressed_genes = sc_up_genes

supp_tab['activated_genes'] = supp_tab['genes'].isin(activated_genes)
supp_tab['repressed_genes'] = supp_tab['genes'].isin(repressed_genes)   

# save the dataframe 
os.makedirs("../output/supplementary_table", exist_ok=True)
supp_tab.to_csv("../output/supplementary_table/sg_creba_network.csv")

