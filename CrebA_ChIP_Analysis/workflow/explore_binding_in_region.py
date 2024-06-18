import pandas as pd
import numpy as np 
import os

###### to make the working directory ###### 
output_path = "../output/match_manual_automated_in_region_peaks"
if os.path.isdir(output_path) == False: 
    os.makedirs(output_path)    

##### this is the function to reproduce data ##### 
def find_bound_in_region_genes(peak_path, DE_path):
    narrow_peaks = pd.read_csv(peak_path, index_col=0)
    narrow_peaks = narrow_peaks.loc[narrow_peaks['in_region_gene'].isna() == False, :]
    bound_genes = list()
    for i in narrow_peaks.index: 
        bound_genes = bound_genes + narrow_peaks.loc[i, 'in_region_gene'].split(",")
        
    output_df = pd.DataFrame(index = ['manual_bound_all', 'manual_unbound_all', 'manual_bound_SG', 'manual_unbound_SG'], columns = ['down_genes', 'up_genes', 'static_genes'])
    for i in range(0, 3): 
        if i == 0: 
            type_gene = 'down_genes'
        elif i == 1: 
            type_gene = 'up_genes'
        elif i == 2: 
            type_gene = 'static_genes'

        DE_df = pd.read_excel(DE_path, sheet_name = i)
        DE_df['bound_status'] = 0
        DE_df.loc[DE_df['target_genes'].isin(bound_genes), 'bound_status'] = 1

        sub_DE = DE_df.loc[[True if 'Unbound' in x else False for x in DE_df['Bound Status']], :]
        output_df.loc['manual_unbound_all', type_gene] = np.sum(sub_DE['bound_status']) / sub_DE.shape[0]

        sub_DE = sub_DE.loc[sub_DE['Salivary Gland'] != 'other', :]
        output_df.loc['manual_unbound_SG', type_gene] = np.sum(sub_DE['bound_status']) / sub_DE.shape[0]

        sub_DE = DE_df.loc[[True if 'Bound' in x else False for x in DE_df['Bound Status']], :]
        output_df.loc['manual_bound_all', type_gene] = np.sum(sub_DE['bound_status']) / sub_DE.shape[0]

        sub_DE = sub_DE.loc[sub_DE['Salivary Gland'] != 'other', :]
        output_df.loc['manual_bound_SG', type_gene] = np.sum(sub_DE['bound_status']) / sub_DE.shape[0]

    return output_df 

spcg_df = pd.read_excel('../input/SPCG_files/SPCG List.xlsx')

###### +/- 1kb fkh sage intersect genes ######
peak_path = "../output/match_nearest_gene/fkh_sage_intersect_500_genes.csv"
DE_path = "../input/manual_curated_DE_genes/manual_curated_DE.xlsx"
output_df = find_bound_in_region_genes(peak_path, DE_path)
print(output_df)

narrow_peaks = pd.read_csv(peak_path, index_col=0)
narrow_peaks = narrow_peaks.loc[narrow_peaks['in_region_gene'].isna() == False, :]
print(narrow_peaks.shape[0])
bound_genes = list()
for i in narrow_peaks.index: 
    bound_genes = bound_genes + narrow_peaks.loc[i, 'in_region_fly_id'].split(",")
bound_genes = np.unique(bound_genes)
print(len(bound_genes))

spcg_df['fkh_sage_intersect_peaks_bound'] = False
spcg_df.loc[spcg_df['Drosophila FBgn'].isin(bound_genes), 'fkh_sage_intersect_peaks_bound'] = True
np.sum(spcg_df['Drosophila FBgn'].isin(bound_genes)) / spcg_df.shape[0]

##### +/- 1kb fkh sage unique ##### 
peak_path = "../output/match_nearest_gene/fkh_sage_unique_500_genes.csv"
DE_path = "../input/manual_curated_DE_genes/manual_curated_DE.xlsx"
output_df = find_bound_in_region_genes(peak_path, DE_path)
print(output_df)

narrow_peaks = pd.read_csv(peak_path, index_col=0)
narrow_peaks = narrow_peaks.loc[narrow_peaks['in_region_gene'].isna() == False, :]
print(narrow_peaks.shape[0])
bound_genes = list()
for i in narrow_peaks.index: 
    bound_genes = bound_genes + narrow_peaks.loc[i, 'in_region_fly_id'].split(",")
bound_genes = np.unique(bound_genes)
print(len(bound_genes))

spcg_df['fkh_sage_unique_peaks_bound'] = False
spcg_df.loc[spcg_df['Drosophila FBgn'].isin(bound_genes), 'fkh_sage_unique_peaks_bound'] = True
np.sum(spcg_df['Drosophila FBgn'].isin(bound_genes)) / spcg_df.shape[0]

spcg_df.to_csv(os.path.join(output_path, "spcg_match.csv"))