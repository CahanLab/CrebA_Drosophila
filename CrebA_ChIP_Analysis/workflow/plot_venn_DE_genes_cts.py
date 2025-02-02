import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3_circles
import pandas as pd 
import numpy as np 
import os 
import itertools
import venn

# make the color palette 
color_palette = dict()
color_palette['Salivary Gland'] = '#1f77b4'
color_palette['Amnioserosa'] = '#ff7f0e'
color_palette['Plasmatocytes'] = '#2ca02c'
color_palette['Fat Body'] = '#9467bd'
color_palette['Trachea'] = '#d62728'  

# make the dictionary of DE genes 
DE_genes_dict = dict()
DE_genes_up_dict = dict()
DE_genes_down_dict = dict()

for tmp_key in color_palette.keys(): 
    DE_genes_dict[tmp_key] = pd.read_csv("../../analysis/results/v19/DE_genes_early_crebA_wt/" + tmp_key + "/mut_DE_genes.csv", index_col = 0)
    DE_genes_dict[tmp_key] = DE_genes_dict[tmp_key].loc[DE_genes_dict[tmp_key]['pval'] < 0.05, :].copy()
    DE_genes_dict[tmp_key] = DE_genes_dict[tmp_key].loc[np.abs(DE_genes_dict[tmp_key]['logFC']) > 0.15, :].copy()
    DE_genes_down_dict[tmp_key] = set(DE_genes_dict[tmp_key].loc[DE_genes_dict[tmp_key]['logFC'] < 0, 'feature'])
    DE_genes_up_dict[tmp_key] = set(DE_genes_dict[tmp_key].loc[DE_genes_dict[tmp_key]['logFC'] > 0, 'feature'])


output_path = "../output/ct_genes_comparison"
if os.path.isdir(output_path) == False:
    os.makedirs(output_path)

ct_combinations = list(itertools.combinations(color_palette.keys(), 2))

# get the down regulated genes 
for tmp_combo in ct_combinations: 
    # get the down regulated genes 
    ax = venn2([set(DE_genes_down_dict[tmp_combo[0]]), set(DE_genes_down_dict[tmp_combo[1]])], 
      (tmp_combo[0], tmp_combo[1]), 
      set_colors = [color_palette[tmp_combo[0]], color_palette[tmp_combo[1]]], alpha = 0.8)  
    for patch in ax.patches:
        if patch:  # Only modify existing patches (some could be None)
            patch.set_edgecolor('black')   # Set border color
            patch.set_linewidth(1.5)         # Set border width
    plt.title(tmp_combo[0] + " and " + tmp_combo[1] + " down genes comparisons") 
    plt.savefig(os.path.join(output_path, tmp_combo[0] + "_" + tmp_combo[1] + '_down.png'), dpi=300, bbox_inches='tight')
    plt.clf()

# get the up regulated genes 
for tmp_combo in ct_combinations: 
    # get the down regulated genes 
    ax = venn2([set(DE_genes_up_dict[tmp_combo[0]]), set(DE_genes_up_dict[tmp_combo[1]])], 
      (tmp_combo[0], tmp_combo[1]), 
      set_colors = [color_palette[tmp_combo[0]], color_palette[tmp_combo[1]]], alpha = 0.8)  
    for patch in ax.patches:
        if patch:  # Only modify existing patches (some could be None)
            patch.set_edgecolor('black')   # Set border color
            patch.set_linewidth(1.5)         # Set border width
    plt.title(tmp_combo[0] + " and " + tmp_combo[1] + " up genes comparisons") 
    plt.savefig(os.path.join(output_path, tmp_combo[0] + "_" + tmp_combo[1] + '_up.png'), dpi=300, bbox_inches='tight')
    plt.clf()

# make the 5 overlaps 
color_palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd', '#d62728']
venn.venn(DE_genes_down_dict, cmap = color_palette, fontsize=13, legend_loc = 'upper right')
plt.title("Down genes overlap across 5 tissues")
plt.savefig(os.path.join(output_path, 'all_5cts_down.png'), bbox_inches='tight', dpi=600)
plt.clf()

venn.venn(DE_genes_up_dict, cmap = color_palette, fontsize=13)
plt.title("Up genes overlap across 5 tissues")
plt.savefig(os.path.join(output_path, 'all_5cts_up.png'), dpi=600, bbox_inches='tight')
plt.clf()