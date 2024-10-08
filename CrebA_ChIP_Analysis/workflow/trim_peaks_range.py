import pandas as pd
import numpy as np 
import os
import argparse 

parser = argparse.ArgumentParser(description='trim the range of peaks')
parser.add_argument('--range', type=int, help = 'the length of the trim', default = 250)
parser.add_argument('--bed', type=str, help = 'the path for the bed files')
parser.add_argument('--out', type=str, help = 'the output path')
args = parser.parse_args()

reg_range = args.range 
input_bed = args.bed
output_bed = args.out

##### look at all the peaks that are unique enough ###### 
peak_file = pd.read_csv(input_bed, sep = '\t', header=None)
new_start = np.array(peak_file[1] + peak_file[9] - reg_range)
new_start[new_start < 0] = 0
new_end = peak_file[1] + peak_file[9] + reg_range
peak_file[1] = new_start
peak_file[2] = new_end
peak_file[9] = new_end - reg_range - new_start 
peak_file.to_csv(output_bed, sep='\t', header=None, index=False)
