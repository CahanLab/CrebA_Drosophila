import pandas as pd 
import numpy as np 
import os 
import argparse 

parser = argparse.ArgumentParser(description='match genes to nearest gene')
parser.add_argument('--bed', type=str, help = 'the path for the bed files')
parser.add_argument('--maxLength', type = int, help = 'maximum distance for matching nearest gene')
parser.add_argument('--out', type=str, help = 'the output path')
args = parser.parse_args()

bed_input = args.bed 
bed_output = args.out
max_length = args.maxLength

##### load in the required data 
# look at the intersecting peaks 500bp limited 
narrow_peaks = pd.read_csv(os.path.join(bed_input), sep="\t", header=None)
promoter_tss = pd.read_csv("../output/tss_table/tss_table.txt")

###### create the nearest distance 
narrow_peaks['nearest_fly_id_1'] = None 
narrow_peaks['nearest_gene_1'] = None 
narrow_peaks['nearest_distance_1'] = None

narrow_peaks['nearest_fly_id_2'] = None 
narrow_peaks['nearest_gene_2'] = None 
narrow_peaks['nearest_distance_2'] = None

narrow_peaks['nearest_fly_id'] = None
narrow_peaks['nearest_gene'] = None

narrow_peaks[0] = narrow_peaks[0].str.removeprefix('chr')

def find_smallest_positive(lst, max_length = 2000):
    positive_numbers = [num for num in lst if num > 0]    
    if positive_numbers:  # Check if there are any positive numbers
        smallest_positive = np.min(positive_numbers)
        if smallest_positive > max_length:
            return None
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
    
    #sub_promoter_tss['down_dist'] = (narrow_peaks.loc[i, 2] - 250) - sub_promoter_tss['tss']
    #sub_promoter_tss['up_dist'] = sub_promoter_tss['tss'] - (narrow_peaks.loc[i, 2] - 250)

    sub_promoter_tss['down_dist'] = (int(narrow_peaks.loc[i, 1]) + int(narrow_peaks.loc[i, 9])) - sub_promoter_tss['tss']
    sub_promoter_tss['up_dist'] = sub_promoter_tss['tss'] - (int(narrow_peaks.loc[i, 1]) + int(narrow_peaks.loc[i, 9]))
    
    smallest_pos_dist = find_smallest_positive(list(sub_promoter_tss['down_dist']), max_length)
    if smallest_pos_dist != None: 
        fly_id = sub_promoter_tss.loc[sub_promoter_tss['down_dist'] == smallest_pos_dist, 'gene_id'].str.cat(sep = ",")
        gene_id = sub_promoter_tss.loc[sub_promoter_tss['down_dist'] == smallest_pos_dist, 'gene_name'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'nearest_fly_id_1'] = fly_id
        narrow_peaks.loc[i, 'nearest_gene_1'] = gene_id
        narrow_peaks.loc[i, 'nearest_distance_1'] = smallest_pos_dist

    smallest_pos_dist = find_smallest_positive(list(sub_promoter_tss['up_dist']), max_length)
    if smallest_pos_dist != None: 
        fly_id = sub_promoter_tss.loc[sub_promoter_tss['up_dist'] == smallest_pos_dist, 'gene_id'].str.cat(sep = ",")
        gene_id = sub_promoter_tss.loc[sub_promoter_tss['up_dist'] == smallest_pos_dist, 'gene_name'].str.cat(sep = ",")
        narrow_peaks.loc[i, 'nearest_fly_id_2'] = fly_id
        narrow_peaks.loc[i, 'nearest_gene_2'] = gene_id
        narrow_peaks.loc[i, 'nearest_distance_2'] = smallest_pos_dist
    
narrow_peaks.to_csv(bed_output)

