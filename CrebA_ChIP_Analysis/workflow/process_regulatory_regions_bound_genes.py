import pandas as pd
import numpy as np 
import os 

output_path = '../output/bound_genes_FIMO'
os.makedirs(output_path, exist_ok=True)

# load in the data 
tss_tab = pd.read_csv("../output/tss_table/tss_table.txt")
range_size = 1500 

# load in the genes that are bound genes 
down_bound_genes_df = pd.read_csv("../output/find_bound_DE_genes/down_DE.csv", index_col = 0)
up_bound_genes_df = pd.read_csv("../output/find_bound_DE_genes/up_DE.csv", index_col = 0)

down_bound_genes_df = down_bound_genes_df.loc[down_bound_genes_df['bound'] == True, :].copy()
down_bound_genes_df = down_bound_genes_df.loc[down_bound_genes_df['SC_DE'] == True, :].copy()

up_bound_genes_df = up_bound_genes_df.loc[up_bound_genes_df['bound'] == True, :].copy()
up_bound_genes_df = up_bound_genes_df.loc[up_bound_genes_df['SC_DE'] == True, :].copy()

# make the down bed files 
sub_tss_tab = tss_tab.loc[tss_tab['gene_name'].isin(down_bound_genes_df['genes']), :].copy()
down_bed = pd.DataFrame()
down_bed[0] = sub_tss_tab['chrom']
down_bed[1] = sub_tss_tab['tss'] - range_size
down_bed.loc[down_bed[1] < 1, 1] = 1
