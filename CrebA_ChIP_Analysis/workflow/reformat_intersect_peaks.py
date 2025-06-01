import pandas as pd
import numpy as np 
import os
import argparse 

parser = argparse.ArgumentParser(description='reformat the peaks')
parser.add_argument('--bed', type=str, help = 'the path for the bed files')
parser.add_argument('--out', type=str, help = 'the output path')
args = parser.parse_args()

input_bed = args.bed
output_bed = args.out

##### look at all the peaks that are unique enough ###### 
peak_file = pd.read_csv(input_bed, sep = '\t', header=None)
#peak_file[3] = np.array(['intersect_peaks_' + str(x) for x in peak_file.index])
new_mid = (peak_file[2] - peak_file[1]) / 2
peak_file[9] = new_mid.astype(int)
peak_file.to_csv(output_bed, sep='\t', header=None, index=False)
