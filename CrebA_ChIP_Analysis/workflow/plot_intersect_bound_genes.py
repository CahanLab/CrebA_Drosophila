import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
import pandas as pd 
import numpy as np 
import os 

# create the directory if needed 
output_path = "../output/bound_spcg_DE_venn"
if os.path.isdir(output_path) == False:
    os.makedirs(output_path)

down_DE_genes_df = pd.read_csv("../output/find_bound_DE_genes/down_DE.csv", index_col=0)
down_genes = np.unique(down_DE_genes_df['genes'])
down_DE_genes_df = down_DE_genes_df.loc[down_DE_genes_df['SC_DE'] == True, :]
sc_down_genes = np.unique(down_DE_genes_df['genes'])

up_DE_genes_df = pd.read_csv("../output/find_bound_DE_genes/up_DE.csv", index_col=0)
up_genes = np.unique(up_DE_genes_df['genes'])
up_DE_genes_df = up_DE_genes_df.loc[up_DE_genes_df['SC_DE'] == True, :]
sc_up_genes = np.unique(up_DE_genes_df['genes'])

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

###### there are to make plots #######
venn3([set(bound_genes), set(spcgs_genes), set(down_genes)], ('Bound Genes', 'SPCGs', 'Down Genes'))  
plt.title("Bound Genes, SPCGs, Down Genes") 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_down_venn.png'), dpi=300, bbox_inches='tight')

venn3([set(bound_genes), set(spcgs_genes), set(sc_down_genes)], ('Bound Genes', 'SPCGs', 'SC Down Genes'))  
plt.title("Bound Genes, SPCGs, SC Down Genes") 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_scdown_venn.png'), dpi=300, bbox_inches='tight')

venn3([set(bound_genes), set(spcgs_genes), set(up_genes)], ('Bound Genes', 'SPCGs', 'Up Genes'))  
plt.title("Bound Genes, SPCGs, Up Genes") 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_up_venn.png'), dpi=300, bbox_inches='tight')

venn3([set(bound_genes), set(spcgs_genes), set(sc_up_genes)], ('Bound Genes', 'SPCGs', 'SC Up Genes'))  
plt.title("Bound Genes, SPCGs, SC Up Genes") 
plt.savefig(os.path.join(output_path, 'bound_SPCGs_scup_venn.png'), dpi=300, bbox_inches='tight')

venn3([set(bound_genes), set(down_genes), set(up_genes)], ('Bound Genes', 'Down Genes', 'Up Genes'))  
plt.title("Bound Genes, Down Genes, Up Genes") 
plt.savefig(os.path.join(output_path, 'bound_down_up_venn.png'), dpi=300, bbox_inches='tight')

