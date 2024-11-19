import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import pandas as pd 
import numpy as np 
import os 

# create the directory if needed 
output_path = "../output/bound_spcg_DE_venn_4cts"
if os.path.isdir(output_path) == False:
    os.makedirs(output_path)

flybase_df = pd.read_csv("../input/flybase_gene_conversion/conversion_tab.csv")
flybase_df.index = flybase_df['flybase']

down_DE_genes_df = pd.read_csv("../output/find_bound_DE_genes_4cts/down_DE.csv", index_col=0)
down_genes = np.unique(down_DE_genes_df['genes'])

sc_down_DE_genes_df = down_DE_genes_df.loc[down_DE_genes_df['SC_DE'] == True, :]
sc_down_genes = np.unique(sc_down_DE_genes_df['genes'])

ma_down_DE_genes_df = down_DE_genes_df.loc[down_DE_genes_df['MA_DE'] == True, :]
MA_down_genes = np.unique(ma_down_DE_genes_df['genes'])

insitu_down_DE_genes_df = down_DE_genes_df.loc[down_DE_genes_df['in_situ_DE'] == True, :]
insitu_genes = np.unique(insitu_down_DE_genes_df['genes'])

# get the genes that are up 
up_DE_genes_df = pd.read_csv("../output/find_bound_DE_genes_4cts/up_DE.csv", index_col=0)
up_genes = np.unique(up_DE_genes_df['genes'])

sc_up_DE_genes_df = up_DE_genes_df.loc[up_DE_genes_df['SC_DE'] == True, :]
sc_up_genes = np.unique(sc_up_DE_genes_df['genes'])

ma_up_DE_genes_df = up_DE_genes_df.loc[up_DE_genes_df['MA_DE'] == True, :]
MA_up_genes = np.unique(ma_up_DE_genes_df['genes'])

# get the spcgs 
spcgs_df = pd.read_csv("../output/find_bound_SPCGs/spcg_match_1000.csv", index_col=0)
spcgs_genes = np.unique(spcgs_df['new_gene_id'])

# get all the bound genes 
crebA_bound_df = pd.read_csv("../output/match_nearest_gene/oregon_intersect_genes_1000.csv", index_col = 0)
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

###### there are to make plots #######
# ma, sc, insitus
ax = venn3([set(MA_down_genes), set(sc_down_genes), set(insitu_genes)], 
      ('MA down genes', 'SC down genes', 'in-situ'), 
      set_colors = ['#6B717E', '#66c2a5', '#8DA0CB'], alpha = 0.8)  
plt.title("MA down, SC down, SPCG") 
plt.savefig(os.path.join(output_path, 'ma_down-sc_down-insitu.png'), dpi=300, bbox_inches='tight')
plt.clf()
# comment -- seems like in-situ and sc down genes seem to be the best for just SG 

venn3([set(bound_genes), set(spcgs_genes), set(down_genes)], 
      ('Bound genes', 'SPCGs', 'Activated genes'), 
      set_colors = ['#e78ac3', '#ffd92f', '#66c2a5'], alpha = 0.8)  
plt.title("CrebA bound and activated genes") 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_down_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

venn3([set(bound_genes), set(spcgs_genes), set(sc_down_genes)], 
      ('Bound genes', 'SPCGs', 'Activated genes'), 
      set_colors = ['#e78ac3', '#ffd92f', '#66c2a5'], alpha = 0.8)  
plt.title("CrebA bound and activated genes in 4 tissues") 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_down_sc_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

venn3([set(bound_genes), set(spcgs_genes), set(np.concatenate((sc_down_genes, insitu_genes)))], 
      ('Bound genes', 'SPCGs', 'Activated genes'), 
      set_colors = ['#e78ac3', '#ffd92f', '#66c2a5'], alpha = 0.8)  
plt.title("CrebA bound and activated genes in 4 tissues") 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_down_sc_insitu_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

venn3([set(bound_genes), set(spcgs_genes), set(up_genes)], 
      ('Bound genes', 'SPCGs', 'Repressed Genes'), 
      set_colors = ['#e78ac3', '#ffd92f', '#fc8d62'], alpha = 0.8)  
plt.title("CrebA bound and repressed genes") 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_up_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

venn3([set(bound_genes), set(spcgs_genes), set(sc_up_genes)], 
      ('Bound genes', 'SPCGs', 'Repressed Genes'), 
      set_colors = ['#e78ac3', '#ffd92f', '#fc8d62'], alpha = 0.8)  
plt.title("CrebA bound and repressed genes in 4 tissues") 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_up_sc_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

venn3([set(bound_genes), set(down_genes), set(up_genes)], ('Bound Genes', 'Down Genes', 'Up Genes'))  
plt.title("Bound Genes, Down Genes, Up Genes") 
plt.savefig(os.path.join(output_path, 'bound_down_up_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

# ma, sc, spcg
fig, ax = plt.subplots(figsize=(8, 8))  # width and height in inches
ax = venn3([set(MA_down_genes), set(sc_down_genes), set(spcgs_genes)], 
      ('MA down genes', 'SC down genes', 'SPCGs'), 
      set_colors = ['#6B717E', '#66c2a5', '#ffd92f'], alpha = 0.8)  
plt.title("MA down, SC down, SPCG") 
plt.savefig(os.path.join(output_path, 'ma_down-sc_down-SPCGs.png'), dpi=300, bbox_inches='tight')
plt.clf()

# ma up, sc up spcg
fig, ax = plt.subplots(figsize=(8, 8))  # width and height in inches
ax = venn3([set(MA_up_genes), set(sc_up_genes), set(spcgs_genes)], 
      ('MA up genes', 'SC up genes', 'SPCGs'), 
      set_colors = ['#B8D1A2', '#fc8d62', '#ffd92f'], alpha = 0.8)  
plt.title("MA up, SC up, SPCG") 
plt.savefig(os.path.join(output_path, 'ma_up-sc_up-SPCGs.png'), dpi=300, bbox_inches='tight')
plt.clf()

