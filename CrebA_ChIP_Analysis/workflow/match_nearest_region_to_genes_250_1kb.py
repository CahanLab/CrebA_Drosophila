import pandas as pd 
import numpy as np 
import os 

# look at the intersecting peaks 500bp limited 
narrow_peaks = pd.read_csv(os.path.join("../output/500_peak_regions/oregon_fkh_sage_intersect_500.narrowPeak"), sep="\t", header=None)
promoter_tss = pd.read_csv("../output/tss_table/tss_table.txt")

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
    start_pos = int(narrow_peaks.loc[i, 2]) - 205 - 1000
    end_pos = int(narrow_peaks.loc[i, 2]) - 250 + 1000

    sub_promoter_tss = sub_promoter_tss.loc[sub_promoter_tss['chrom'] == temp_chrom, :].copy()
    sub_promoter_tss['in_region'] = np.logical_and(sub_promoter_tss['tss'] >= start_pos, sub_promoter_tss['tss'] <= end_pos)
    
    if sub_promoter_tss.shape[0] == 0: 
        next
    
    if np.sum(sub_promoter_tss['in_region']) > 0: 
        narrow_peaks.loc[i, 'in_region_fly_id'] = sub_promoter_tss.loc[sub_promoter_tss['in_region'] == True, 'gene_id'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'in_region_gene'] = sub_promoter_tss.loc[sub_promoter_tss['in_region'] == True, 'gene_name'].str.cat(sep = ",")
    
narrow_peaks.to_csv("../output/match_nearest_gene/oregon_fkh_sage_intersect_in_region_genes.csv")

##### load in the unique peaks ######
narrow_peaks = pd.read_csv(os.path.join("../output/500_peak_regions/oregon_fkh_sage_unique_500.narrowPeak"), sep="\t", header=None)
promoter_tss = pd.read_csv("../output/tss_table/tss_table.txt")

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
    start_pos = int(narrow_peaks.loc[i, 2]) - 205 - 1000
    end_pos = int(narrow_peaks.loc[i, 2]) - 250 + 1000

    sub_promoter_tss = sub_promoter_tss.loc[sub_promoter_tss['chrom'] == temp_chrom, :].copy()
    sub_promoter_tss['in_region'] = np.logical_and(sub_promoter_tss['tss'] >= start_pos, sub_promoter_tss['tss'] <= end_pos)
    
    if sub_promoter_tss.shape[0] == 0: 
        next
    
    if np.sum(sub_promoter_tss['in_region']) > 0: 
        narrow_peaks.loc[i, 'in_region_fly_id'] = sub_promoter_tss.loc[sub_promoter_tss['in_region'] == True, 'gene_id'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'in_region_gene'] = sub_promoter_tss.loc[sub_promoter_tss['in_region'] == True, 'gene_name'].str.cat(sep = ",")
    
narrow_peaks.to_csv("../output/match_nearest_gene/oregon_fkh_sage_unique_in_region_genes.csv")

##### load in the unique peaks ######
narrow_peaks = pd.read_csv(os.path.join("../output/500_peak_regions/CrebA_oregon_500.narrowPeak"), sep="\t", header=None)
promoter_tss = pd.read_csv("../output/tss_table/tss_table.txt")

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
    start_pos = int(narrow_peaks.loc[i, 2]) - 205 - 1000
    end_pos = int(narrow_peaks.loc[i, 2]) - 250 + 1000

    sub_promoter_tss = sub_promoter_tss.loc[sub_promoter_tss['chrom'] == temp_chrom, :].copy()
    sub_promoter_tss['in_region'] = np.logical_and(sub_promoter_tss['tss'] >= start_pos, sub_promoter_tss['tss'] <= end_pos)
    
    if sub_promoter_tss.shape[0] == 0: 
        next
    
    if np.sum(sub_promoter_tss['in_region']) > 0: 
        narrow_peaks.loc[i, 'in_region_fly_id'] = sub_promoter_tss.loc[sub_promoter_tss['in_region'] == True, 'gene_id'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'in_region_gene'] = sub_promoter_tss.loc[sub_promoter_tss['in_region'] == True, 'gene_name'].str.cat(sep = ",")
    
narrow_peaks.to_csv("../output/match_nearest_gene/oregon_in_region_genes.csv")

##### load in the fkh_sage unique peaks ######
narrow_peaks = pd.read_csv(os.path.join("../output/500_peak_regions/fkh_sage_unique_500.narrowPeak"), sep="\t", header=None)
promoter_tss = pd.read_csv("../output/tss_table/tss_table.txt")

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
    start_pos = int(narrow_peaks.loc[i, 2]) - 205 - 1000
    end_pos = int(narrow_peaks.loc[i, 2]) - 250 + 1000

    sub_promoter_tss = sub_promoter_tss.loc[sub_promoter_tss['chrom'] == temp_chrom, :].copy()
    sub_promoter_tss['in_region'] = np.logical_and(sub_promoter_tss['tss'] >= start_pos, sub_promoter_tss['tss'] <= end_pos)
    
    if sub_promoter_tss.shape[0] == 0: 
        next
    
    if np.sum(sub_promoter_tss['in_region']) > 0: 
        narrow_peaks.loc[i, 'in_region_fly_id'] = sub_promoter_tss.loc[sub_promoter_tss['in_region'] == True, 'gene_id'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'in_region_gene'] = sub_promoter_tss.loc[sub_promoter_tss['in_region'] == True, 'gene_name'].str.cat(sep = ",")
    
narrow_peaks.to_csv("../output/match_nearest_gene/fkh_sage_unique_in_region_genes.csv")

