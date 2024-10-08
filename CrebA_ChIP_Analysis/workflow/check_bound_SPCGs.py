import pandas as pd
import numpy as np 

spcg_1000 = pd.read_csv("../output/find_bound_SPCGs/spcg_match_1000.csv")
spcg_1000 = spcg_1000.loc[spcg_1000['fkh_sage_intersect_peaks_bound'] == False, :]

spcg_2000 = pd.read_csv("../output/find_bound_SPCGs/spcg_match_2000.csv")
spcg_2000 = spcg_2000.loc[spcg_2000['fkh_sage_intersect_peaks_bound'] == False, :]

np.setdiff1d(spcg_1000['Drosophila Gene'], spcg_2000['Drosophila Gene'])