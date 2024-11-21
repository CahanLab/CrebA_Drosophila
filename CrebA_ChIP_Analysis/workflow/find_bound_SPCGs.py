import pandas as pd
import numpy as np 
import os
import argparse 

parser = argparse.ArgumentParser(description='find bound SPCGs')
parser.add_argument('--input', type=str, help = 'input file')
parser.add_argument('--output', type=str, help = 'output file')
args = parser.parse_args()

spcg_df = pd.read_excel('../input/SPCG_files/SPCG List.xlsx')
conversion_tab = pd.read_csv("../input/flybase_gene_conversion/conversion_tab.csv")
conversion_tab.index = conversion_tab['flybase']
sub_conversion_tab = conversion_tab.loc[spcg_df['Drosophila FBgn'], :]
spcg_df['new_gene_id'] = np.array(sub_conversion_tab['gene_names'])

###### +/- 1kb fkh sage intersect genes ######
peak_path = args.input
output_path = args.output

narrow_peaks = pd.read_csv(peak_path, index_col=0)
concatenated_array = np.concatenate((np.array(narrow_peaks['nearest_fly_id_1'].dropna()), np.array(narrow_peaks['nearest_fly_id_2'].dropna())))
concatenated_array = np.unique(concatenated_array)
concatenated_array = list(concatenated_array)
for element in concatenated_array:
    # Split the string by whitespace
    split_elements = element.split(",")
    if len(split_elements) > 1:
        concatenated_array.extend(split_elements)
        concatenated_array.remove(element)

bound_genes = concatenated_array

spcg_df['peaks_bound'] = False
spcg_df.loc[spcg_df['Drosophila FBgn'].isin(bound_genes), 'peaks_bound'] = True
spcg_df.to_csv(output_path, index=False)