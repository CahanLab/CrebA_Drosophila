import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import pandas as pd 
import numpy as np 
import os 

# create the directory if needed 
output_path = "../output/bound_spcg_DE_venn"
if os.path.isdir(output_path) == False:
    os.makedirs(output_path)

flybase_df = pd.read_csv("../input/flybase_gene_conversion/conversion_tab.csv")
flybase_df.index = flybase_df['flybase']

down_DE_genes_df = pd.read_csv("../output/find_bound_DE_genes/down_DE.csv", index_col=0)
down_genes = np.unique(down_DE_genes_df['genes'])
down_DE_genes_df = down_DE_genes_df.loc[down_DE_genes_df['SC_DE'] == True, :]
sc_down_genes = np.unique(down_DE_genes_df['genes'])

up_DE_genes_df = pd.read_csv("../output/find_bound_DE_genes/up_DE.csv", index_col=0)
up_genes = np.unique(up_DE_genes_df['genes'])
up_DE_genes_df = up_DE_genes_df.loc[up_DE_genes_df['SC_DE'] == True, :]
sc_up_genes = np.unique(up_DE_genes_df['genes'])

MA_down_genes_df = pd.read_csv("../input/nate_ma_DE/down_ma_DE_NateTranslated.csv")
MA_down_genes = MA_down_genes_df['Converted Gene Name']

MA_up_genes_df = pd.read_csv("../input/nate_ma_DE/up_ma_DE_TranslatedNate.csv")
MA_up_genes = MA_up_genes_df['Translated Gene Names']

# get the spcgs 
spcgs_df = pd.read_csv("../output/find_bound_SPCGs/spcg_match.csv", index_col=0)
spcgs_genes = np.unique(spcgs_df['new_gene_id'])

# get all the bound genes 
crebA_bound_df = pd.read_csv("../output/match_nearest_gene/fkh_sage_intersect_genes.csv", index_col = 0)
crebA_bound_df = crebA_bound_df.loc[crebA_bound_df['in_region_gene'].isna() == False, :]
bound_genes = list()
for i in crebA_bound_df.index: 
    bound_genes = bound_genes + crebA_bound_df.loc[i, 'in_region_gene'].split(",")
bound_genes = np.unique(bound_genes)

# see how many of the genes are uniquely bounded 
unique_bound_genes = list()
for my_index in crebA_bound_df.index: 
    in_region_genes = crebA_bound_df.loc[my_index, 'in_region_gene']
    if "," not in in_region_genes: 
        unique_bound_genes.append(in_region_genes)
multi_bound_genes = np.setdiff1d(bound_genes, unique_bound_genes)

###### there are to make plots #######
venn3([set(bound_genes), set(spcgs_genes), set(down_genes)], 
      ('Bound genes', 'SPCGs', 'Activated genes'), 
      set_colors = ['#e78ac3', '#ffd92f', '#66c2a5'], alpha = 0.8)  
plt.title("CrebA bound and activated genes") 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_down_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

venn3([set(bound_genes), set(spcgs_genes), set(sc_down_genes)], 
      ('Bound genes', 'SPCGs', 'Activated genes'), 
      set_colors = ['#e78ac3', '#ffd92f', '#66c2a5'], alpha = 0.8)  
plt.title("CrebA bound and activated genes in salivary gland") 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_down_venn_sg.png'), dpi=300, bbox_inches='tight')
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
plt.title("CrebA bound and repressed genes in salivary gland") 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_up_venn_sg.png'), dpi=300, bbox_inches='tight')
plt.clf()

venn3([set(bound_genes), set(down_genes), set(up_genes)], ('Bound Genes', 'Down Genes', 'Up Genes'))  
plt.title("Bound Genes, Down Genes, Up Genes") 
plt.savefig(os.path.join(output_path, 'bound_down_up_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

# multi single activate genes in sg 
fig, ax = plt.subplots(figsize=(10, 10))  # width and height in inches
ax = venn3([set(multi_bound_genes), set(unique_bound_genes), set(sc_down_genes)], 
      ('Multi-bound genes', 'Single-bound genes', 'Activated genes'), 
      set_colors = ['#FF99CC', '#C060A6', '#66c2a5'], alpha = 0.8)  
plt.title("mutli, single activated genes in SG") 
plt.savefig(os.path.join(output_path, 'multi-single-activate-sg.png'), dpi=300, bbox_inches='tight')
plt.clf()

# multi single repressed genes in sg 
fig, ax = plt.subplots(figsize=(10, 10))  # width and height in inches
ax = venn3([set(multi_bound_genes), set(unique_bound_genes), set(sc_up_genes)], 
      ('Multi-bound genes', 'Single-bound genes', 'Repressed genes'), 
      set_colors = ['#FF99CC', '#C060A6', '#fc8d62'], alpha = 0.8)  
plt.title("mutli, single repressed genes in SG") 
plt.savefig(os.path.join(output_path, 'multi-single-repressed-sg.png'), dpi=300, bbox_inches='tight')
plt.clf()

# multi single in SPCGs
fig, ax = plt.subplots(figsize=(10, 10))  # width and height in inches
ax = venn3([set(multi_bound_genes), set(unique_bound_genes), set(spcgs_genes)], 
      ('Multi-bound genes', 'Single-bound genes', 'SPCGs'), 
      set_colors = ['#FF99CC', '#C060A6', '#ffd92f'], alpha = 0.8)  
plt.title("mutli, single bound in SPCGs") 
plt.savefig(os.path.join(output_path, 'multi-single-SPCGs.png'), dpi=300, bbox_inches='tight')
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