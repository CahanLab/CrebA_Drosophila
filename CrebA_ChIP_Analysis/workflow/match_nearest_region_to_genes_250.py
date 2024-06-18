import pandas as pd 
import numpy as np 
import os 

# look at the intersecting peaks 500bp limited 
narrow_peaks = pd.read_csv(os.path.join("../output/500_peak_regions/oregon_fkh_sage_intersect_500.narrowPeak"), sep="\t", header=None)
promoter_tss = pd.read_csv("../output/tss_table/tss_table.txt")

narrow_peaks['nearest_fly_id_1'] = None 
narrow_peaks['nearest_gene_1'] = None 
narrow_peaks['nearest_distance_1'] = None

narrow_peaks['nearest_fly_id_2'] = None 
narrow_peaks['nearest_gene_2'] = None 
narrow_peaks['nearest_distance_2'] = None

narrow_peaks['nearest_fly_id'] = None
narrow_peaks['nearest_gene'] = None

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
    

    sub_promoter_tss['down_dist'] = (narrow_peaks.loc[i, 2] - 250) - sub_promoter_tss['tss']
    sub_promoter_tss['up_dist'] = sub_promoter_tss['tss'] - (narrow_peaks.loc[i, 2] - 250)
    sub_promoter_tss['abs_dist'] = np.abs(sub_promoter_tss['down_dist'])
    smallest_pos_dist = find_smallest_positive(list(sub_promoter_tss['down_dist']))

    if smallest_pos_dist != None: 
        fly_id = sub_promoter_tss.loc[sub_promoter_tss['down_dist'] == smallest_pos_dist, 'gene_id'].str.cat(sep = ",")
        gene_id = sub_promoter_tss.loc[sub_promoter_tss['down_dist'] == smallest_pos_dist, 'gene_name'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'nearest_fly_id_1'] = fly_id
        narrow_peaks.loc[i, 'nearest_gene_1'] = gene_id
        narrow_peaks.loc[i, 'nearest_distance_1'] = smallest_pos_dist

    smallest_pos_dist = find_smallest_positive(list(sub_promoter_tss['up_dist']))
    if smallest_pos_dist != None: 
        fly_id = sub_promoter_tss.loc[sub_promoter_tss['up_dist'] == smallest_pos_dist, 'gene_id'].str.cat(sep = ",")
        gene_id = sub_promoter_tss.loc[sub_promoter_tss['up_dist'] == smallest_pos_dist, 'gene_name'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'nearest_fly_id_2'] = fly_id
        narrow_peaks.loc[i, 'nearest_gene_2'] = gene_id
        narrow_peaks.loc[i, 'nearest_distance_2'] = smallest_pos_dist

    smallest_pos_dist = find_smallest_positive(list(sub_promoter_tss['abs_dist']))
    if smallest_pos_dist != None:
        fly_id = sub_promoter_tss.loc[sub_promoter_tss['abs_dist'] == smallest_pos_dist, 'gene_id'].str.cat(sep = ",")
        gene_id = sub_promoter_tss.loc[sub_promoter_tss['abs_dist'] == smallest_pos_dist, 'gene_name'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'nearest_fly_id'] = fly_id
        narrow_peaks.loc[i, 'nearest_gene'] = gene_id
    
if os.path.isdir("../output/match_nearest_gene/") == False:
    os.makedirs("../output/match_nearest_gene/")
narrow_peaks.to_csv("../output/match_nearest_gene/oregon_fkh_sage_intersect_nearest_genes.csv")

##### load in the unique peaks ######
narrow_peaks = pd.read_csv(os.path.join("../output/500_peak_regions/oregon_fkh_sage_unique_500.narrowPeak"), sep="\t", header=None)
promoter_tss = pd.read_csv("../output/tss_table/tss_table.txt")

narrow_peaks['nearest_fly_id_1'] = None 
narrow_peaks['nearest_gene_1'] = None 
narrow_peaks['nearest_distance_1'] = None

narrow_peaks['nearest_fly_id_2'] = None 
narrow_peaks['nearest_gene_2'] = None 
narrow_peaks['nearest_distance_2'] = None

narrow_peaks['nearest_fly_id'] = None
narrow_peaks['nearest_gene'] = None

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
    

    sub_promoter_tss['down_dist'] = (narrow_peaks.loc[i, 2] - 250) - sub_promoter_tss['tss']
    sub_promoter_tss['up_dist'] = sub_promoter_tss['tss'] - (narrow_peaks.loc[i, 2] - 250)
    sub_promoter_tss['abs_dist'] = np.abs(sub_promoter_tss['down_dist'])
    smallest_pos_dist = find_smallest_positive(list(sub_promoter_tss['down_dist']))

    if smallest_pos_dist != None: 
        fly_id = sub_promoter_tss.loc[sub_promoter_tss['down_dist'] == smallest_pos_dist, 'gene_id'].str.cat(sep = ",")
        gene_id = sub_promoter_tss.loc[sub_promoter_tss['down_dist'] == smallest_pos_dist, 'gene_name'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'nearest_fly_id_1'] = fly_id
        narrow_peaks.loc[i, 'nearest_gene_1'] = gene_id
        narrow_peaks.loc[i, 'nearest_distance_1'] = smallest_pos_dist

    smallest_pos_dist = find_smallest_positive(list(sub_promoter_tss['up_dist']))
    if smallest_pos_dist != None: 
        fly_id = sub_promoter_tss.loc[sub_promoter_tss['up_dist'] == smallest_pos_dist, 'gene_id'].str.cat(sep = ",")
        gene_id = sub_promoter_tss.loc[sub_promoter_tss['up_dist'] == smallest_pos_dist, 'gene_name'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'nearest_fly_id_2'] = fly_id
        narrow_peaks.loc[i, 'nearest_gene_2'] = gene_id
        narrow_peaks.loc[i, 'nearest_distance_2'] = smallest_pos_dist

    smallest_pos_dist = find_smallest_positive(list(sub_promoter_tss['abs_dist']))
    if smallest_pos_dist != None:
        fly_id = sub_promoter_tss.loc[sub_promoter_tss['abs_dist'] == smallest_pos_dist, 'gene_id'].str.cat(sep = ",")
        gene_id = sub_promoter_tss.loc[sub_promoter_tss['abs_dist'] == smallest_pos_dist, 'gene_name'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'nearest_fly_id'] = fly_id
        narrow_peaks.loc[i, 'nearest_gene'] = gene_id
    
if os.path.isdir("../output/match_nearest_gene/") == False:
    os.makedirs("../output/match_nearest_gene/")
narrow_peaks.to_csv("../output/match_nearest_gene/oregon_fkh_sage_unique_nearest_genes.csv")

##### load in the unique peaks ######
narrow_peaks = pd.read_csv(os.path.join("../output/500_peak_regions/CrebA_oregon_500.narrowPeak"), sep="\t", header=None)
promoter_tss = pd.read_csv("../output/tss_table/tss_table.txt")

narrow_peaks['nearest_fly_id_1'] = None 
narrow_peaks['nearest_gene_1'] = None 
narrow_peaks['nearest_distance_1'] = None

narrow_peaks['nearest_fly_id_2'] = None 
narrow_peaks['nearest_gene_2'] = None 
narrow_peaks['nearest_distance_2'] = None

narrow_peaks['nearest_fly_id'] = None
narrow_peaks['nearest_gene'] = None

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
    

    sub_promoter_tss['down_dist'] = (narrow_peaks.loc[i, 2] - 250) - sub_promoter_tss['tss']
    sub_promoter_tss['up_dist'] = sub_promoter_tss['tss'] - (narrow_peaks.loc[i, 2] - 250)
    sub_promoter_tss['abs_dist'] = np.abs(sub_promoter_tss['down_dist'])
    smallest_pos_dist = find_smallest_positive(list(sub_promoter_tss['down_dist']))

    if smallest_pos_dist != None: 
        fly_id = sub_promoter_tss.loc[sub_promoter_tss['down_dist'] == smallest_pos_dist, 'gene_id'].str.cat(sep = ",")
        gene_id = sub_promoter_tss.loc[sub_promoter_tss['down_dist'] == smallest_pos_dist, 'gene_name'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'nearest_fly_id_1'] = fly_id
        narrow_peaks.loc[i, 'nearest_gene_1'] = gene_id
        narrow_peaks.loc[i, 'nearest_distance_1'] = smallest_pos_dist

    smallest_pos_dist = find_smallest_positive(list(sub_promoter_tss['up_dist']))
    if smallest_pos_dist != None: 
        fly_id = sub_promoter_tss.loc[sub_promoter_tss['up_dist'] == smallest_pos_dist, 'gene_id'].str.cat(sep = ",")
        gene_id = sub_promoter_tss.loc[sub_promoter_tss['up_dist'] == smallest_pos_dist, 'gene_name'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'nearest_fly_id_2'] = fly_id
        narrow_peaks.loc[i, 'nearest_gene_2'] = gene_id
        narrow_peaks.loc[i, 'nearest_distance_2'] = smallest_pos_dist

    smallest_pos_dist = find_smallest_positive(list(sub_promoter_tss['abs_dist']))
    if smallest_pos_dist != None:
        fly_id = sub_promoter_tss.loc[sub_promoter_tss['abs_dist'] == smallest_pos_dist, 'gene_id'].str.cat(sep = ",")
        gene_id = sub_promoter_tss.loc[sub_promoter_tss['abs_dist'] == smallest_pos_dist, 'gene_name'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'nearest_fly_id'] = fly_id
        narrow_peaks.loc[i, 'nearest_gene'] = gene_id
    
if os.path.isdir("../output/match_nearest_gene/") == False:
    os.makedirs("../output/match_nearest_gene/")
narrow_peaks.to_csv("../output/match_nearest_gene/oregon_nearest_genes.csv")

