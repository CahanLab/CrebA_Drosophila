import pandas as pd 
import numpy as np 
import os 

# create the directory if needed 
if os.path.isdir("../output/match_nearest_gene/") == False:
    os.makedirs('../output/match_nearest_gene/')

promoter_tss = pd.read_csv("../output/tss_table/tss_table.txt")

# match the peaks to genes 
for tmp_file in ['fkh_sage_intersect', 'fkh_sage_unique']:
    narrow_peaks = pd.read_csv(os.path.join("../output/range_peak_regions/" + tmp_file + ".narrowPeak"), sep="\t", header=None)
    narrow_peaks['in_region_fly_id'] = None 
    narrow_peaks['in_region_gene'] = None 
    narrow_peaks[0] = narrow_peaks[0].str.removeprefix('chr')

    def find_smallest_positive(lst):
        positive_numbers = [num for num in lst if num > 0]    
        if positive_numbers:  # Check if there are any positive numbers
            smallest_positive = np.min(positive_numbers)
            return smallest_positive
        else:
            return None
        
    for i in narrow_peaks.index:
        sub_promoter_tss = promoter_tss.copy()
        temp_chrom = narrow_peaks.loc[i, 0]
        start_pos = int(narrow_peaks.loc[i, 1]) 
        end_pos = int(narrow_peaks.loc[i, 2]) 
        sub_promoter_tss = sub_promoter_tss.loc[sub_promoter_tss['chrom'] == temp_chrom, :].copy()
        sub_promoter_tss['in_region'] = np.logical_and(sub_promoter_tss['tss'] >= start_pos, sub_promoter_tss['tss'] <= end_pos)
        if sub_promoter_tss.shape[0] == 0: 
            next
        if np.sum(sub_promoter_tss['in_region']) > 0: 
            narrow_peaks.loc[i, 'in_region_fly_id'] = sub_promoter_tss.loc[sub_promoter_tss['in_region'] == True, 'gene_id'].str.cat(sep = ",")
            narrow_peaks.loc[i, 'in_region_gene'] = sub_promoter_tss.loc[sub_promoter_tss['in_region'] == True, 'gene_name'].str.cat(sep = ",")
    narrow_peaks.to_csv("../output/match_nearest_gene/" + tmp_file + "_genes.csv")

