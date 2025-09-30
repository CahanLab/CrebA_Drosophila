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
parser.add_argument('--MA_include', type=str, help = 'whether to include MA', default='No')
parser.add_argument('--MA_counter_file', type=str, help = 'whether to include MA', default='')

args = parser.parse_args()

bed_input = args.bed 
DE_genes = args.DE_genes
spcg_input = args.spcg
output_path = args.out 
MA_include = args.MA_include
MA_counter_file = args.MA_counter_file

# this is the main script 
bed_file = pd.read_csv(bed_input, index_col=0)
bed_file = bed_file.loc[np.logical_and(bed_file['nearest_fly_id_1'].isna(), bed_file['nearest_fly_id_2'].isna()) == False, :]

SG_de_genes = pd.read_csv("../input/DE_genes_early_crebA_wt/Salivary Gland/mut_DE_genes.csv", index_col = 0)

def check_intersect_genes(bed_file, gene_list):
    output_array = []
    for tmp_index in bed_file.index: 
        geneset_1 = str(bed_file.loc[tmp_index, 'nearest_gene_1']).split(",")
        geneset_2 = str(bed_file.loc[tmp_index, 'nearest_gene_2']).split(",")
        geneset = geneset_1 + geneset_2
        if len(np.intersect1d(geneset, gene_list)) == 0: 
            output_array.append(False)
        else: 
            output_array.append(True)
    return np.array(output_array)

for tmp_index in bed_file.index: 
    tmp_name_list = [bed_file.loc[tmp_index, '3'], bed_file.loc[tmp_index, 'nearest_gene_1'], bed_file.loc[tmp_index, 'nearest_gene_2']]
    tmp_name_list = [item for item in tmp_name_list if not isinstance(item, float) or not math.isnan(item)]
    new_name = "-".join(tmp_name_list)
    bed_file.loc[tmp_index, '3'] = new_name

if spcg_input == '' and DE_genes == '':
    trimmed_bed = bed_file.iloc[:, 0:10]
elif spcg_input == '' and DE_genes != "":
    DE_genes_df = pd.read_csv(DE_genes, index_col=0)

    if 'down_DE' in DE_genes:
        if MA_include == 'No':
            DE_genes_df = DE_genes_df.loc[np.logical_or(DE_genes_df['SC_DE'] == True, DE_genes_df['in_situ_DE'] == True), :]
        elif MA_include == 'MA_total':
            counter_DE_genes_df = pd.read_csv(MA_counter_file, index_col=0)
            i_genes = np.intersect1d(counter_DE_genes_df['genes'], DE_genes_df['genes'])
            DE_genes_df = DE_genes_df.loc[DE_genes_df['genes'].isin(i_genes) == False, :]
        elif MA_include == 'MA_SG':
            SG_genes = pd.read_csv("../input/early_sg_DE_genes/markers_genes_early_SG.csv", index_col = 0)
            SG_genes = SG_genes.loc[SG_genes['pct.1'] >= 0.1, :]
            
            mut_diff_genes = pd.read_csv("../input/DE_genes_early_crebA_wt/Salivary Gland/mut_DE_genes.csv", index_col = 0)
            mut_diff_genes = mut_diff_genes.loc[mut_diff_genes['logFC'] < 0, :]

            MA_SG = np.intersect1d(SG_genes.index, mut_diff_genes['feature'])

            # SC or in situ DE 
            SC_insitu_index = np.logical_or(DE_genes_df['SC_DE'] == True, DE_genes_df['in_situ_DE'] == True)
            # MA DE and SG genes
            MA_SG_index = np.logical_and(DE_genes_df['MA_DE'] == True, DE_genes_df['genes'].isin(MA_SG))
            # either first condition or second condition
            all_index = np.logical_or(SC_insitu_index, MA_SG_index)

            DE_genes_df = DE_genes_df.loc[all_index, :]
    else:
        if MA_include != 'MA_total':
            DE_genes_df = DE_genes_df.loc[DE_genes_df['SC_DE'] == True, :]
        else:
            counter_DE_genes_df = pd.read_csv(MA_counter_file, index_col=0)
            i_genes = np.intersect1d(counter_DE_genes_df['genes'], DE_genes_df['genes'])
            DE_genes_df = DE_genes_df.loc[DE_genes_df['genes'].isin(i_genes) == False, :]

    sub_bed_file = bed_file.loc[check_intersect_genes(bed_file, DE_genes_df['genes']), :]
    trimmed_bed = sub_bed_file.iloc[:, 0:10]
elif spcg_input != '' and DE_genes == '':
    spcg_file = pd.read_excel(spcg_input)
    sub_bed_file = bed_file.loc[np.logical_or(bed_file['nearest_fly_id_1'].isin(spcg_file['Drosophila FBgn']), 
                                            bed_file['nearest_fly_id_2'].isin(spcg_file['Drosophila FBgn'])), :]
    trimmed_bed = sub_bed_file.iloc[:, 0:10]

# rank sub_bed_file 
SG_de_genes.index = SG_de_genes['feature']

sub_bed_file['logFC'] = None

for tmp_index in sub_bed_file.index:
    tmp_gene_1 = sub_bed_file.loc[tmp_index, 'nearest_gene_1']
    tmp_gene_2 = sub_bed_file.loc[tmp_index, 'nearest_gene_2']

    if pd.isna(tmp_gene_1): 
        tmp_gene_1 = 'NO_GENE'
    if pd.isna(tmp_gene_2):
        tmp_gene_2 = 'NO_GENE'
    list_genes = tmp_gene_1.split(",") + tmp_gene_2.split(",")
    
    logFC_list = list()
    for tmp_gene in list_genes:
        if tmp_gene in SG_de_genes.index:
            logFC_list.append(SG_de_genes.loc[tmp_gene, 'logFC'])
    if len(logFC_list) == 0:
        # if the gene is only in MA and in-situs but does not exist in SC 
        sub_bed_file.loc[tmp_index, 'logFC'] = 99
    else:
        sub_bed_file.loc[tmp_index, 'logFC'] = np.min(logFC_list)

sub_bed_file = sub_bed_file.sort_values(by='logFC', ascending=True)
sub_bed_file.to_csv(output_path+'.csv', index = False)

trimmed_bed = sub_bed_file.iloc[:, 0:10]
trimmed_bed.to_csv(output_path, header=False, sep='\t', index=False)

