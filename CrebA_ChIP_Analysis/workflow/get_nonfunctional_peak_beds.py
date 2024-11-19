import pandas as pd 
import numpy as np 
import os 
import argparse 
import math

parser = argparse.ArgumentParser(description='match genes to nearest gene')
parser.add_argument('--bed', type=str, help = 'the path for the bed files')
parser.add_argument('--DE_genes_down', type = str, help = 'the path for DE down genes', default='')
parser.add_argument('--DE_genes_up', type = str, help = 'the path for DE up genes', default='')
parser.add_argument('--scRNA_results', type = str, help = 'the path for single cell RNA-seq DE results', default='')
parser.add_argument('--num_genes', type = int, help = 'the number of nonfunctional genes', default=100)

parser.add_argument('--out', type=str, help = 'the output path')
args = parser.parse_args()

bed_input = args.bed 
DE_genes_down = args.DE_genes_down
DE_genes_up = args.DE_genes_up
scRNA_results = args.scRNA_results
num_genes = args.num_genes

output_path = args.out 

# this is the main script 
bed_file = pd.read_csv(bed_input, index_col=0)
bed_file = bed_file.loc[np.logical_and(bed_file['nearest_fly_id_1'].isna(), bed_file['nearest_fly_id_2'].isna()) == False, :]
all_bound_genes = np.unique(list(bed_file['nearest_gene_1']) + list(bed_file['nearest_gene_2']))

concatenated_array = np.concatenate((np.array(bed_file['nearest_gene_1'].dropna()), np.array(bed_file['nearest_gene_2'].dropna())))
concatenated_array = np.unique(concatenated_array)
concatenated_array = list(concatenated_array)
for element in concatenated_array:
    # Split the string by whitespace
    split_elements = element.split(",")
    if len(split_elements) > 1:
        concatenated_array.extend(split_elements)
        concatenated_array.remove(element)
bound_genes = np.unique(concatenated_array)

for tmp_index in bed_file.index: 
    tmp_name_list = [bed_file.loc[tmp_index, '3'], bed_file.loc[tmp_index, 'nearest_gene_1'], bed_file.loc[tmp_index, 'nearest_gene_2']]
    tmp_name_list = [item for item in tmp_name_list if not isinstance(item, float) or not math.isnan(item)]
    new_name = "-".join(tmp_name_list)
    bed_file.loc[tmp_index, '3'] = new_name

DE_down_df = pd.read_csv(DE_genes_down, index_col=0)
DE_down_df = DE_down_df.loc[DE_down_df['bound'] == True, :] # all the genes have at least 1 DE modality 

DE_up_df = pd.read_csv(DE_genes_up, index_col = 0)
DE_up_df = DE_up_df.loc[DE_up_df['bound'] == True, :]

functional_genes = np.unique(list(DE_down_df['genes']) + list(DE_up_df['genes']))
if scRNA_results == '':
    sub_bed_file = bed_file.loc[np.logical_and(bed_file['nearest_gene_1'].isin(functional_genes) == False, 
                                                bed_file['nearest_gene_2'].isin(functional_genes) == False), :]
else: 
    scRNA_DE = pd.read_csv(scRNA_results, index_col=0)
    scRNA_DE = scRNA_DE.loc[scRNA_DE.index.isin(functional_genes) == False, :] # not funcitonal genes
    scRNA_DE = scRNA_DE.loc[scRNA_DE.index.isin(bound_genes), :] # intersect with bound genes
    scRNA_DE = scRNA_DE.sort_values(by='avg_log2FC', key=abs) # rank by log2FC 
    sub_scRNA_DE = scRNA_DE.iloc[0:num_genes, :] # get least top DE genes
    sub_bed_file = bed_file.loc[np.logical_or(bed_file['nearest_gene_1'].isin(sub_scRNA_DE.index), 
                                            bed_file['nearest_gene_2'].isin(sub_scRNA_DE.index)), :]

trimmed_bed = sub_bed_file.iloc[:, 0:10]
trimmed_bed.to_csv(output_path, header=False, sep='\t', index=False)