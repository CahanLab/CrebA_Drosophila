import pandas as pd
import numpy as np 
import os 
import random

output_path = '../output/bound_genes_regulatory_regions_400bp_null'
os.makedirs(output_path, exist_ok=True)

##### get the maximum length for each chromosome
def get_max_length(file_path, chr_prefix = False):
    output_dict = dict()
    with open(file_path, 'r') as file: 
        for line in file: 
            if line.startswith(">"):
                if chr_prefix == False:
                    output_dict[line.split(" ")[0][1:]] = line.split(" ")[6].removeprefix('length=').removesuffix(";")
                else: 
                    output_dict['chr' + line.split(" ")[0][1:]] = line.split(" ")[6].removeprefix('length=').removesuffix(";")
    return output_dict
max_dict = get_max_length('../input/reference_genome/dmel-all-chromosome-r6.33.fasta')

# load in the data 
reg_region_tab = pd.read_csv("../output/tss_table/regulatory_regions_400bp.bed", sep = '\t', header=None)
all_bound_genes = pd.read_csv("../output/plot_location_crebA_binding/all_bound_genes_dist.csv", index_col=0)
target_bed = reg_region_tab.loc[reg_region_tab[6].isin(all_bound_genes['bound_gene']) == False, :].copy()
bound_bed = reg_region_tab.loc[reg_region_tab[6].isin(all_bound_genes['bound_gene']) == True, :].copy()

reg_tab_bed = target_bed.loc[:, [0, 2, 3, 6]].copy()
reg_tab_bed = reg_tab_bed.loc[reg_tab_bed[6].isna() == False, :]
reg_tab_bed.to_csv(os.path.join(output_path, 'no_bound_genes_reg_regions.bed'), sep = '\t', header = None, index = None)

##### randomly sample 400 bp regions ###### 
np.random.seed(1234)
random.seed(1234)

num_trials = 2000
seq_size = 400
elements = np.array(bound_bed[0].value_counts().index)
probabilities = np.array(bound_bed[0].value_counts() / np.sum(bound_bed[0].value_counts()))

rand_bed_tab = pd.DataFrame(index=list(range(0, num_trials)), columns = [0, 1, 2, 3])
for i in range(0, num_trials):
    sampled_chrom = np.random.choice(elements, p=probabilities)
    rand_bed_tab.loc[i, 0] = sampled_chrom
    rand_coords = list(range(0, int(max_dict[sampled_chrom])))
    rand_bed_tab.loc[i, 1] = random.choice(rand_coords)
    rand_bed_tab.loc[i, 2] = rand_bed_tab.loc[i, 1] + seq_size
    rand_bed_tab.loc[i, 3] = 'random_sample_' + str(i)
    print(i)
rand_bed_tab.to_csv(os.path.join(output_path, 'random_reg_regions.bed'), sep = '\t', header = None, index = None)
