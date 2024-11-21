import pandas as pd
import numpy as np 
import argparse
import os

parser = argparse.ArgumentParser(description='get all the bound genes genes')
parser.add_argument('--DE_genes', type=str, help = 'folder path for DE genes')
parser.add_argument('--bound_genes', type=str, help = 'the path for the bound genes')
parser.add_argument('--out', type=str, help = 'the output path')
args = parser.parse_args()

DE_genes_path = args.DE_genes
bound_genes_file = args.bound_genes

output_folder = args.out

os.makedirs(output_folder, exist_ok=True)

down_DE_genes_df = pd.read_csv(os.path.join(DE_genes_path, "down_DE.csv"), index_col=0)
up_DE_genes_df = pd.read_csv(os.path.join(DE_genes_path, "up_DE.csv"), index_col=0)

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
concatenated_array = np.unique(concatenated_array)

compiled_df = pd.DataFrame()
compiled_df['bound_genes'] = concatenated_array
compiled_df.index = compiled_df['bound_genes']

# grab all the down genes 
compiled_df['MA_DE_down'] = False
down_genes = np.array(down_DE_genes_df.loc[np.logical_and(down_DE_genes_df['MA_DE'] == True, down_DE_genes_df['bound'] == True), 'genes'])
compiled_df.loc[down_genes, 'MA_DE_down'] = True

compiled_df['SC_DE_down'] = False
down_genes = np.array(down_DE_genes_df.loc[np.logical_and(down_DE_genes_df['SC_DE'] == True, down_DE_genes_df['bound'] == True), 'genes'])
compiled_df.loc[down_genes, 'SC_DE_down'] = True

compiled_df['in_situ_DE_down'] = False
down_genes = np.array(down_DE_genes_df.loc[np.logical_and(down_DE_genes_df['in_situ_DE'] == True, down_DE_genes_df['bound'] == True), 'genes'])
compiled_df.loc[down_genes, 'in_situ_DE_down'] = True

# grab all the up genes 
compiled_df['MA_DE_up'] = False
up_genes = np.array(up_DE_genes_df.loc[np.logical_and(up_DE_genes_df['MA_DE'] == True, up_DE_genes_df['bound'] == True), 'genes'])
compiled_df.loc[up_genes, 'MA_DE_up'] = True

compiled_df['SC_DE_up'] = False
up_genes = np.array(up_DE_genes_df.loc[np.logical_and(up_DE_genes_df['SC_DE'] == True, up_DE_genes_df['bound'] == True), 'genes'])
compiled_df.loc[up_genes, 'SC_DE_up'] = True

compiled_df.to_csv(os.path.join(output_folder, "all_bound_genes_data.csv"), index=False)