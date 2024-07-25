import pandas as pd 
import numpy as np 
import os
import matplotlib.pyplot as plt
import scipy.stats as stats 
from statsmodels.stats.multitest import fdrcorrection

##### get the fimo results for up and down ###### 
input_path = '../output/bound_genes_regulatory_regions'
down_TF = pd.read_csv(os.path.join(input_path, 'down_flyfactorsurvey_hits.csv'), index_col = 0)
up_TF = pd.read_csv(os.path.join(input_path, 'up_flyfactorsurvey_hits.csv'), index_col = 0)

##### analyze the data ######
down_TF = down_TF.loc[up_TF.index, :]

TF_motif_list = []
odd_ratio_list = []
pval_list = []
for tmp_TF in down_TF.index: 
    pos_up = np.sum(up_TF.loc[tmp_TF, :] == 1)
    neg_up = np.sum(up_TF.loc[tmp_TF, :] == 0)
    pos_down = np.sum(down_TF.loc[tmp_TF, :] == 1)
    neg_down = np.sum(down_TF.loc[tmp_TF, :] == 0)
    cont_tab = [[pos_up, neg_up], [pos_down, neg_down]]
    odd_ratio, p_value = stats.fisher_exact(cont_tab) 
    TF_motif_list.append(tmp_TF)
    odd_ratio_list.append(odd_ratio)
    pval_list.append(p_value)

output_df = pd.DataFrame()
output_df['TF_motif'] = TF_motif_list
output_df['odd_ratio'] = odd_ratio_list
output_df['p_value'] = pval_list
output_df = output_df.sort_values(by = 'p_value')
output_df['fdr_val'] = fdrcorrection(np.array(output_df['p_value']))[1]