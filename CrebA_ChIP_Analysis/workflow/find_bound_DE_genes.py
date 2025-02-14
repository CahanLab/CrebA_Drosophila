import pandas as pd
import numpy as np 
import datetime
import os
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

import argparse 

parser = argparse.ArgumentParser(description='find bound DE genes')
parser.add_argument('--boundGenes', type=str, help = 'the path for bound genes')
parser.add_argument(
    '--cts', 
    nargs='+',  # '+' means one or more arguments
    type=str,   # Specify the type for each item in the list
    help='A list of celltypes'
)
parser.add_argument('--include_late_sg', type=str, help = 'whether to include the late SG genes', default = 'no')
parser.add_argument('--out', type=str, help = 'the output path')
args = parser.parse_args()

bound_genes_file = args.boundGenes 
cts_list = args.cts
output_path = args.out
include_late_sg = args.include_late_sg

# create the directory if needed 
if os.path.isdir(output_path) == False:
    os.makedirs(output_path)

# make the dictionary of DE genes 
DE_genes_dict = dict()
DE_genes_up_dict = dict()
DE_genes_down_dict = dict()

for tmp_key in cts_list: 
    DE_genes_dict[tmp_key] = pd.read_csv("../../analysis/results/v19/DE_genes_early_crebA_wt/" + tmp_key + "/mut_DE_genes.csv", index_col = 0)
    DE_genes_dict[tmp_key] = DE_genes_dict[tmp_key].loc[DE_genes_dict[tmp_key]['pval'] < 0.05, :].copy()
    DE_genes_dict[tmp_key] = DE_genes_dict[tmp_key].loc[np.abs(DE_genes_dict[tmp_key]['logFC']) > 0.15, :].copy()
    DE_genes_down_dict[tmp_key] = np.array(DE_genes_dict[tmp_key].loc[DE_genes_dict[tmp_key]['logFC'] < 0, 'feature'])
    DE_genes_up_dict[tmp_key] = np.array(DE_genes_dict[tmp_key].loc[DE_genes_dict[tmp_key]['logFC'] > 0, 'feature'])

if include_late_sg == 'yes':
    tmp_key = 'late_Salivary Gland'
    cts_list.append(tmp_key)
    DE_genes_dict[tmp_key] = pd.read_csv("../../analysis/results/v19/DE_genes_crebA_wt/" + "Salivary Gland/mut_DE_genes.csv", index_col = 0)
    DE_genes_dict[tmp_key] = DE_genes_dict[tmp_key].loc[DE_genes_dict[tmp_key]['pval'] < 0.05, :].copy()
    DE_genes_dict[tmp_key] = DE_genes_dict[tmp_key].loc[np.abs(DE_genes_dict[tmp_key]['logFC']) > 0.15, :].copy()
    DE_genes_down_dict[tmp_key] = np.array(DE_genes_dict[tmp_key].loc[DE_genes_dict[tmp_key]['logFC'] < 0, 'feature'])
    DE_genes_up_dict[tmp_key] = np.array(DE_genes_dict[tmp_key].loc[DE_genes_dict[tmp_key]['logFC'] > 0, 'feature'])

# get all down genes, up genes and the unique ones 
down_genes = list() 
up_genes = list() 
for tmp_key in cts_list: 
    down_genes = down_genes + list(DE_genes_down_dict[tmp_key])
    up_genes = up_genes + list(DE_genes_up_dict[tmp_key])

# plot out the intersect between all up and down genes 
ax = venn2([set(down_genes), set(up_genes)], 
    ('Down Genes', 'Up Genes'), 
    set_colors = ["#7B596B", "#4A697F"], alpha = 0.8)  
for patch in ax.patches:
    if patch:  # Only modify existing patches (some could be None)
        patch.set_edgecolor('black')   # Set border color
        patch.set_linewidth(1.5)         # Set border width
plt.savefig(os.path.join(output_path, 'up_down_comparison.png'), dpi=300, bbox_inches='tight')

# look at what genes among the 44 genes are conflicting 
conflict_genes = np.intersect1d(down_genes, up_genes)

conflict_df = pd.DataFrame(columns=['gene', 'down_ct', 'up_ct'])
for tmp_gene in conflict_genes: 
    up_cts = list()
    down_cts = list()
    for tmp_ct in cts_list: 
        if tmp_gene in DE_genes_down_dict[tmp_ct]: 
            down_cts.append(tmp_ct)
        if tmp_gene in DE_genes_up_dict[tmp_ct]:
            up_cts.append(tmp_ct)
    conflict_df.loc[len(conflict_df)] = [tmp_gene, ','.join(down_cts), ','.join(up_cts)]
conflict_df.to_csv(os.path.join(output_path, 'conflict_genes.csv'))

# get all the DE genes from scRNA-seq, microarra and in-situ data 
up_MA_nate = pd.read_excel("../input/nate_ma_DE/CrebA_UpDownDE_GeneList Translated.xlsx", sheet_name = 'UpReg')
down_MA_nate = pd.read_excel("../input/nate_ma_DE/CrebA_UpDownDE_GeneList Translated.xlsx", sheet_name = 'DownReg')

# add in SPCG in-situ data 
flybase_converter = pd.read_csv("../input/flybase_gene_conversion/conversion_tab.csv")
in_situs = pd.read_excel("../input/SPCG_files/SPCG List.xlsx")
in_situs = in_situs.loc[in_situs['require_CrebA_inSitu'] == True, :]
in_situs_genes = np.array(flybase_converter.loc[flybase_converter['flybase'].isin(in_situs['Drosophila FBgn']), 'gene_names'])
reg_genes = dict()

def isNaN(num):
    return num != num

for tmp_type in ['up', 'down']:
    if tmp_type == 'up':
        sub_ma_DE = up_MA_nate
        sc_genes = np.setdiff1d(up_genes, conflict_genes)
    else:
        sub_ma_DE = down_MA_nate
        sc_genes = np.setdiff1d(down_genes, conflict_genes)
    ma_genes = []

    for tmp_gene in list(sub_ma_DE['Gene Symbol']):
        if tmp_gene == datetime.datetime(2009, 9, 5, 0, 0): 
            ma_genes = ma_genes + ['Sep5']
        elif tmp_gene == datetime.datetime(2009, 9, 1, 0, 0): 
            ma_genes = ma_genes + ['Sep1']
        elif tmp_gene == datetime.datetime(2009, 9, 2, 0, 0): 
            ma_genes = ma_genes + ['Sep2']
        elif tmp_gene == datetime.datetime(2009, 12, 1, 0, 0): 
            ma_genes = ma_genes + ['dec-1']
        elif isNaN(tmp_gene) == True: 
            continue
        else:
            ma_genes = ma_genes + [x.strip() for x in tmp_gene.split("///")]
    ma_genes = np.array(ma_genes)
    ma_genes = ma_genes[ma_genes != 'None']
    ma_genes = ma_genes[ma_genes != "---"]
    my_df = pd.DataFrame()
    
    if tmp_type == 'up':
        my_df['genes'] = np.unique(np.concatenate((sc_genes, ma_genes), axis=None))
        my_df['MA_DE'] = my_df['genes'].isin(ma_genes)
        my_df['SC_DE'] = my_df['genes'].isin(sc_genes)
    else: 
        my_df['genes'] = np.unique(np.concatenate((sc_genes, ma_genes, in_situs_genes), axis=None))
        my_df['MA_DE'] = my_df['genes'].isin(ma_genes)
        my_df['SC_DE'] = my_df['genes'].isin(sc_genes)
        my_df['in_situ_DE'] = my_df['genes'].isin(in_situs_genes)
    reg_genes[tmp_type] = my_df

# see if these genes are bound or not 
crebA_bound_df = pd.read_csv(bound_genes_file, index_col = 0)
concatenated_array = np.concatenate((np.array(crebA_bound_df['nearest_gene_1'].dropna()), np.array(crebA_bound_df['nearest_gene_2'].dropna())))
concatenated_array = np.unique(concatenated_array)
concatenated_array = list(concatenated_array)
for element in concatenated_array:
    # Split the string by whitespace
    split_elements = element.split(",")
    if len(split_elements) > 1:
        concatenated_array.extend(split_elements)
        concatenated_array.remove(element)

bound_genes = concatenated_array

for tmp_type in ['up', 'down']:
    reg_genes[tmp_type]['bound'] = False
    reg_genes[tmp_type].loc[reg_genes[tmp_type]['genes'].isin(bound_genes), 'bound'] = True
    reg_genes[tmp_type].to_csv(os.path.join(output_path, tmp_type + "_DE.csv"))
