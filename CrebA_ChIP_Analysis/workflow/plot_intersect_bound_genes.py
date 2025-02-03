import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
import pandas as pd 
import numpy as np 
import os 
import argparse 

parser = argparse.ArgumentParser(description='plot various venn diagrams')
parser.add_argument('--DE_genes', type=str, help = 'the output files with the bound DE genes')
parser.add_argument('--bound_genes', type=str, help = 'the path for the bed files')
parser.add_argument('--out', type=str, help = 'the output path')
parser.add_argument('--sc_type', type=str, help = 'the type of venn diagram to plot for the single cell')
args = parser.parse_args()

DE_genes_path = args.DE_genes
bound_genes_path = args.bound_genes
output_path = args.out
sc_type = args.sc_type

color_palettes = {}
# make an if state to set the color palette 
if sc_type == 'salivary gland':
    color_palettes['repression'] = '#E48482'
    color_palettes['activation'] = '#6BAF92'
    color_palettes['bound'] = '#A4A9D1'
    color_palettes['SPCGs'] = '#D9B76E'
else: 
    color_palettes['repression'] = '#84E4D9'
    color_palettes['activation'] = '#AF6B7F'
    color_palettes['bound'] = '#A9D1A4'
    color_palettes['SPCGs'] = '#7656D9'

if sc_type == 'salivary gland':
    bracket_name = '(SG)'
elif sc_type == 'majority 4cts':
    bracket_name = '(≥3/4 cts)'
else:
    bracket_name = '(' + sc_type + ')'

# create the directory if needed 
if os.path.isdir(output_path) == False:
    os.makedirs(output_path)

flybase_df = pd.read_csv("../input/flybase_gene_conversion/conversion_tab.csv")
flybase_df.index = flybase_df['flybase']

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

###### there are to make plots #######
# ma, sc, insitus
ax = venn3([set(MA_down_genes), set(sc_down_genes), set(insitu_genes)], 
      ('Microarray', 'Single-cell ' + bracket_name, 'in-situ'), 
      set_colors = ['#6B717E', '#66c2a5', '#8DA0CB'], alpha = 0.8)  

# Set font sizes for labels
for label in ax.set_labels:
    label.set_fontsize(16)  # Set size for set labels (A, B)

for label in ax.subset_labels:
    if label:  # Check if the label exists
        label.set_fontsize(14)  # Set size for subset labels
plt.savefig(os.path.join(output_path, 'ma_down-sc_down-insitu.png'), dpi=300, bbox_inches='tight')
plt.clf()
# comment -- seems like in-situ and sc down genes seem to be the best for just SG 


venn3([set(bound_genes), set(spcgs_genes), set(sc_down_genes)], 
      ('Bound genes', 'SPCGs', 'Activated genes'), 
      set_colors = [color_palettes['bound'], color_palettes['SPCGs'], color_palettes['activation']], alpha = 0.8)  
plt.title("CrebA bound and activated genes in " + sc_type) 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_down_sc_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

venn3([set(bound_genes), set(spcgs_genes), set(np.concatenate((sc_down_genes, insitu_genes)))], 
      ('Bound genes', 'SPCGs', 'Activated genes'), 
      set_colors = [color_palettes['bound'], color_palettes['SPCGs'], color_palettes['activation']], alpha = 0.8)  
plt.title("CrebA bound and activated genes in " + sc_type) 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_down_sc_insitu_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

venn2([set(bound_genes),  set(up_genes)], 
      ('Bound genes', 'Repressed Genes'), 
      set_colors = [color_palettes['bound'], color_palettes['repression']], alpha = 0.8)  
plt.title("CrebA bound and repressed genes") 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_up_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

venn2([set(bound_genes), set(sc_up_genes)], 
      ('Bound genes', 'Repressed Genes'), 
      set_colors = [color_palettes['bound'], color_palettes['repression']], alpha = 0.8)  
plt.title("CrebA bound and repressed genes in " + sc_type) 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_up_sc_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

venn3([set(bound_genes), set(down_genes), set(up_genes)], ('Bound Genes', 'Down Genes', 'Up Genes'))  
plt.title("Bound Genes, Down Genes, Up Genes") 
plt.savefig(os.path.join(output_path, 'bound_down_up_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

# ma, sc, spcg
ax = venn3([set(MA_down_genes), set(sc_down_genes), set(spcgs_genes)], 
      ('Microarray', 'Single-cell ' + bracket_name, 'SPCGs'), 
      set_colors = ['#6B717E', '#66c2a5', '#ffd92f'], alpha = 0.8)  
for label in ax.set_labels:
    label.set_fontsize(16)  # Set size for set labels (A, B)

for label in ax.subset_labels:
    if label:  # Check if the label exists
        label.set_fontsize(14)  # Set size for subset labels
plt.savefig(os.path.join(output_path, 'ma_down-sc_down-SPCGs.png'), dpi=300, bbox_inches='tight')
plt.clf()

# ma up, sc up spcg
ax = venn3([set(MA_up_genes), set(sc_up_genes), set(spcgs_genes)], 
      ('Microarray', 'Single-cell ' + bracket_name, 'SPCGs'), 
      set_colors = ['#B8D1A2', '#fc8d62', '#ffd92f'], alpha = 0.8)  
for label in ax.set_labels:
    label.set_fontsize(16)  # Set size for set labels (A, B)

for label in ax.subset_labels:
    if label:  # Check if the label exists
        label.set_fontsize(14)  # Set size for subset labels
plt.savefig(os.path.join(output_path, 'ma_up-sc_up-SPCGs.png'), dpi=300, bbox_inches='tight')
plt.clf()

###### this is the additional analysis that we wanted to do 
# here are the additional code to see if we should include MA entirely or MA and SPCGs
MA_scpg = np.intersect1d(MA_down_genes, spcgs_genes)
venn3([set(bound_genes), set(spcgs_genes), set(np.concatenate((sc_down_genes, insitu_genes, MA_scpg)))], 
      ('Bound genes', 'SPCGs', 'Activated genes'), 
      set_colors = [color_palettes['bound'], color_palettes['SPCGs'], color_palettes['activation']], alpha = 0.8)  
plt.title("CrebA bound and activated genes in " + sc_type) 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_down_sc_insitu_MA-SPCG_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

MA_scpg = np.intersect1d(MA_up_genes, spcgs_genes)

if MA_scpg.size > 0:
    concatenated_genes = np.concatenate((sc_up_genes, MA_scpg))
else:
    concatenated_genes = sc_up_genes

venn2([set(bound_genes), set(concatenated_genes)], 
      ('Bound genes', 'Repressed Genes'), 
      set_colors = [color_palettes['bound'], color_palettes['repression']], alpha = 0.8)  
plt.title("CrebA bound and repressed genes in " + sc_type) 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_up_sc_MA-SPCG_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

# if there is a collision, we go with single-cell and in-situ
sub_MA_down_genes = np.setdiff1d(MA_down_genes, sc_up_genes)

venn3([set(bound_genes), set(spcgs_genes), set(np.concatenate((sc_down_genes, insitu_genes, sub_MA_down_genes)))], 
      ('Bound genes', 'SPCGs', 'Activated genes'), 
      set_colors = [color_palettes['bound'], color_palettes['SPCGs'], color_palettes['activation']], alpha = 0.8)  
plt.title("CrebA bound and activated genes in " + sc_type) 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_down_sc_insitu_MA-sub_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()

sub_MA_up_genes = np.setdiff1d(MA_up_genes, sc_down_genes)
sub_MA_up_genes = np.unique(sub_MA_up_genes[sub_MA_up_genes != ''])

venn2([set(bound_genes), set(np.concatenate((sc_up_genes, sub_MA_up_genes)))], 
      ('Bound genes', 'Repressed Genes'), 
      set_colors = [color_palettes['bound'], color_palettes['repression']], alpha = 0.8)  
plt.title("CrebA bound and repressed genes in " + sc_type) 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_up_sc_MA-sub_venn.png'), dpi=300, bbox_inches='tight')
plt.clf()