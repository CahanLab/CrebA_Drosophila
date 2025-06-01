import pandas as pd
import numpy as np 

sg_no_protein = pd.read_csv("output/supplementary_table/sg_creba_network.csv", index_col = 0)
sg_protein = pd.read_csv("output_protein_only/supplementary_table/sg_creba_network.csv", index_col = 0)

sg_no_protein = sg_no_protein.loc[sg_no_protein['bound'] == True, :]
sg_no_protein = sg_no_protein.loc[sg_no_protein['SPCG'] == True, :]

sg_protein = sg_protein.loc[sg_protein['bound'] == True, :]
sg_protein = sg_protein.loc[sg_protein['SPCG'] == True, :]

np.setdiff1d(sg_protein['genes'], sg_no_protein['genes'])

##### look at the differences between our bound SPCGs and Dorothy's paper ##### 
creba_network = pd.read_csv('output/supplementary_table/sg_creba_network.csv', index_col = 0)

spcg_df = pd.read_excel('input/SPCG_files/SPCG List.xlsx')
conversion_tab = pd.read_csv("input/flybase_gene_conversion/conversion_tab.csv")
conversion_tab.index = conversion_tab['flybase']
sub_conversion_tab = conversion_tab.loc[spcg_df['Drosophila FBgn'], :]
spcg_df['new_gene_id'] = np.array(sub_conversion_tab['gene_names'])
spcg_df = spcg_df.loc[spcg_df['Dorothy_Annotation'] == True, :]

creba_network = creba_network.loc[creba_network['SPCG'] == True, :]
creba_network = creba_network.loc[creba_network['bound'] == True, :]

np.setdiff1d(spcg_df['new_gene_id'], creba_network['genes'])
