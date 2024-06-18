import pandas as pd
import numpy as np 

##### look at the intersecting peaks and limit it to 500 ######
my_intersected_peaks = pd.read_csv("../input/GEO_processed/GSM4213092_oregon_CrebA_peaks.narrowPeak", sep = '\t', header=None)
matt_peaks = pd.read_csv("../input/intersect_peaks/Set1.NEW_OR_fkh_sage_intersection_fromDA.bed", sep = '\t', header = None)
my_intersected_peaks = my_intersected_peaks.loc[my_intersected_peaks[3].isin(matt_peaks[3]), :].copy()

new_start = np.array(my_intersected_peaks[1] + my_intersected_peaks[9] - 250)
new_start[new_start < 0] = 0
new_end = my_intersected_peaks[1] + my_intersected_peaks[9] + 250
my_intersected_peaks[1] = new_start
my_intersected_peaks[2] = new_end

my_intersected_peaks.to_csv("../output/500_peak_regions/oregon_fkh_sage_intersect_500.narrowPeak", sep='\t', header=None, index=False)

##### look at all the peaks that are unique enough ###### 
peak_file = pd.read_csv("../input/GEO_processed/GSM4213092_oregon_CrebA_peaks.narrowPeak", sep = '\t', header=None)
new_start = np.array(peak_file[1] + peak_file[9] - 250)
new_start[new_start < 0] = 0
new_end = peak_file[1] + peak_file[9] + 250
peak_file[1] = new_start
peak_file[2] = new_end
peak_file.to_csv("../output/500_peak_regions/CrebA_oregon_500.narrowPeak", sep='\t', header=None, index=False)

##### look at all the peaks that are unique enough ###### 
peak_file = pd.read_csv("../input/GEO_processed/GSM4213094_fkh_CrebA_peaks.narrowPeak", sep = '\t', header=None)
new_start = np.array(peak_file[1] + peak_file[9] - 250)
new_start[new_start < 0] = 0
new_end = peak_file[1] + peak_file[9] + 250
peak_file[1] = new_start
peak_file[2] = new_end
peak_file.to_csv("../output/500_peak_regions/CrebA_fkh_500.narrowPeak", sep='\t', header=None, index=False)

##### look at all the peaks that are unique enough ###### 
peak_file = pd.read_csv("../input/GEO_processed/GSM4213096_sage_CrebA_peaks.narrowPeak", sep = '\t', header=None)
new_start = np.array(peak_file[1] + peak_file[9] - 250)
new_start[new_start < 0] = 0
new_end = peak_file[1] + peak_file[9] + 250
peak_file[1] = new_start
peak_file[2] = new_end
peak_file.to_csv("../output/500_peak_regions/CrebA_sage_500.narrowPeak", sep='\t', header=None, index=False)