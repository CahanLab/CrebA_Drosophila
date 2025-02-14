import numpy as np 
import pandas as pd
import os

os.makedirs("../output/ranked_genes_fasta", exist_ok=True)

# load in logFC 
sc_DE = pd.read_csv("../../analysis/results/v19/DE_genes_early_crebA_wt/" + 'Salivary Gland' + "/mut_DE_genes.csv", index_col = 0)
sc_DE.index = sc_DE['feature']

# load in TSS 
tss_tab = pd.read_csv("../output/tss_table/tss_table.txt")
tss_tab.index = tss_tab['gene_name']

# read in the DE genes 
DE_genes_df = pd.read_csv("../output/find_bound_DE_genes/down_DE.csv", index_col=0)
DE_genes_df = DE_genes_df.loc[DE_genes_df['bound'] == True, :]
DE_genes_df = DE_genes_df.loc[DE_genes_df['MA_DE'] == True, :]
DE_genes_df = DE_genes_df.loc[DE_genes_df['SC_DE'] == True, :]
DE_genes_df['logFC'] = np.array(sc_DE.loc[DE_genes_df['genes'], 'logFC'])
DE_genes_df = DE_genes_df.sort_values(by='logFC', ascending=True)
DE_genes_df.index = DE_genes_df['genes']
DE_genes_df.to_csv("../output/ranked_genes_fasta/down_genes_rank.csv", sep = '\t', index = False)

bed_df = pd.DataFrame(columns = ['chrom', 'start', 'end', 'gene'])
for tmp_gene in DE_genes_df.index: 
    bed_df.loc[len(bed_df)] = [tss_tab.loc[tmp_gene, 'chrom'], 
                    tss_tab.loc[tmp_gene, 'tss'] - 250, 
                    tss_tab.loc[tmp_gene, 'tss'] + 250, 
                    tmp_gene]
bed_df.to_csv("../output/ranked_genes_fasta/down_genes.bed", sep = '\t', index = False, header = False)

# get the down regulated genes 
DE_genes_df = pd.read_csv("../output/find_bound_DE_genes/up_DE.csv", index_col=0)
DE_genes_df = DE_genes_df.loc[DE_genes_df['bound'] == True, :]
DE_genes_df = DE_genes_df.loc[DE_genes_df['MA_DE'] == True, :]
DE_genes_df = DE_genes_df.loc[DE_genes_df['SC_DE'] == True, :]
DE_genes_df['logFC'] = np.array(sc_DE.loc[DE_genes_df['genes'], 'logFC'])
DE_genes_df = DE_genes_df.sort_values(by='logFC', ascending=False)
DE_genes_df.index = DE_genes_df['genes']
DE_genes_df.to_csv("../output/ranked_genes_fasta/up_genes_rank.csv", sep = '\t', index = False)

bed_df = pd.DataFrame(columns = ['chrom', 'start', 'end', 'gene'])
for tmp_gene in DE_genes_df.index: 
    bed_df.loc[len(bed_df)] = [tss_tab.loc[tmp_gene, 'chrom'], 
                    tss_tab.loc[tmp_gene, 'tss'] - 250, 
                    tss_tab.loc[tmp_gene, 'tss'] + 250, 
                    tmp_gene]
bed_df.to_csv("../output/ranked_genes_fasta/up_genes.bed", sep = '\t', index = False, header = False)

# load top 150 genes 
DE_down_df = pd.read_csv("../output/find_bound_DE_genes/down_DE.csv", index_col=0)
DE_up_df = pd.read_csv("../output/find_bound_DE_genes/up_DE.csv", index_col=0)
functional_genes = np.unique(list(DE_down_df['genes']) + list(DE_up_df['genes']))

bed_file = pd.read_csv('../output/match_nearest_gene/fkh_sage_intersect_genes_1000.csv', index_col=0)
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

scRNA_DE = pd.read_csv("../../analysis/results/v19/DE_genes_early_crebA_wt/" + 'Salivary Gland' + "/mut_DE_genes.csv", index_col = 0)
scRNA_DE.index = scRNA_DE['feature']
scRNA_DE = scRNA_DE.loc[scRNA_DE.index.isin(functional_genes) == False, :] # not funcitonal genes
scRNA_DE = scRNA_DE.loc[scRNA_DE.index.isin(bound_genes), :] # intersect with bound genes
scRNA_DE = scRNA_DE.sort_values(by='logFC', key=abs) # rank by log2FC 
DE_genes_df = scRNA_DE.iloc[0:150, :] # get least top DE genes
DE_genes_df.to_csv("../output/ranked_genes_fasta/non_functional_genes_rank.csv", sep = '\t', index = False)

bed_df = pd.DataFrame(columns = ['chrom', 'start', 'end', 'gene'])
for tmp_gene in DE_genes_df.index: 
    bed_df.loc[len(bed_df)] = [tss_tab.loc[tmp_gene, 'chrom'], 
                    tss_tab.loc[tmp_gene, 'tss'] - 250, 
                    tss_tab.loc[tmp_gene, 'tss'] + 250, 
                    tmp_gene]
bed_df.to_csv("../output/ranked_genes_fasta/150_nonfunctional_genes.bed", sep = '\t', index = False, header = False)
