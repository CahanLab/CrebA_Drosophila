import pandas as pd 
import numpy as np 
import os 
import argparse 
import math

parser = argparse.ArgumentParser(description='match genes to nearest gene')
parser.add_argument('--bed', type=str, help = 'the path for the bed files')
parser.add_argument('--DE_genes', type = str, help = 'the path for DE genes', default='')
parser.add_argument('--spcg', type = str, help = 'the path for SPCGs', default='')
parser.add_argument('--out', type=str, help = 'the output path')
args = parser.parse_args()

bed_input = args.bed 
DE_genes = args.DE_genes
spcg_input = args.spcg
output_path = args.out 

# this is the main script 
bed_file = pd.read_csv(bed_input, index_col=0)
bed_file = bed_file.loc[np.logical_and(bed_file['nearest_fly_id_1'].isna(), bed_file['nearest_fly_id_2'].isna()) == False, :]

for tmp_index in bed_file.index: 
    tmp_name_list = [bed_file.loc[tmp_index, '3'], bed_file.loc[tmp_index, 'nearest_gene_1'], bed_file.loc[tmp_index, 'nearest_gene_2']]
    tmp_name_list = [item for item in tmp_name_list if not isinstance(item, float) or not math.isnan(item)]
    new_name = "-".join(tmp_name_list)
    bed_file.loc[tmp_index, '3'] = new_name

if spcg_input == '' and DE_genes == '':
    trimmed_bed = bed_file.iloc[:, 0:10]
elif spcg_input == '' and DE_genes != "":
    DE_genes_df = pd.read_csv(DE_genes, index_col=0)
    DE_genes_df = DE_genes_df.loc[DE_genes_df['SC_DE'] == True, :]
    sub_bed_file = bed_file.loc[np.logical_or(bed_file['nearest_gene_1'].isin(DE_genes_df['genes']), 
                                              bed_file['nearest_gene_2'].isin(DE_genes_df['genes'])), :]
    trimmed_bed = sub_bed_file.iloc[:, 0:10]
elif spcg_input != '' and DE_genes == '':
    spcg_file = pd.read_excel(spcg_input)
    sub_bed_file = bed_file.loc[np.logical_or(bed_file['nearest_fly_id_1'].isin(spcg_file['Drosophila FBgn']), 
                                            bed_file['nearest_fly_id_2'].isin(spcg_file['Drosophila FBgn'])), :]
    trimmed_bed = sub_bed_file.iloc[:, 0:10]
trimmed_bed.to_csv(output_path, header=False, sep='\t', index=False)