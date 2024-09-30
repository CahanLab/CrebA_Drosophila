import numpy as np 
import pandas as pd
import os

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

gtf_table = pd.read_csv("../input/reference_genome/dmel-all-r6.33.gtf", sep = '\t', header = None)
gtf_table = gtf_table.loc[gtf_table[2] == 'gene', :]

# read in the fasta file to get the length 

max_dict = get_max_length('../input/reference_genome/dmel-all-chromosome-r6.33.fasta')

intersect_peaks = pd.read_csv("../output/range_peak_regions/fkh_sage_intersect.narrowPeak", sep = '\t', header=None)
sage_peaks = pd.read_csv("../input/GEO_processed/GSM4213096_sage_CrebA_peaks.narrowPeak", sep = '\t', header = None)

sage_peaks_new = sage_peaks.loc[sage_peaks[3].isin(intersect_peaks[3]), :]
sage_peaks_new[0] = sage_peaks_new[0].str.removeprefix('chr')

os.makedirs('../output/500bp_regions_intersect_peaks', exist_ok = True)
sage_peaks_new.to_csv("../output/500bp_regions_intersect_peaks/raw_intersect_peaks.bed", sep = '\t', index=None, header=None)