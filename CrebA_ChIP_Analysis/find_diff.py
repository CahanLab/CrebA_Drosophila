import pandas as pd
import numpy as np 

sg_no_protein = pd.read_csv("output/supplementary_table/sg_creba_network.csv", index_col = 0)
sg_protein = pd.read_csv("output_protein_only/supplementary_table/sg_creba_network.csv", index_col = 0)

sg_no_protein = sg_no_protein.loc[sg_no_protein['bound'] == True, :]
sg_no_protein = sg_no_protein.loc[sg_no_protein['SPCG'] == True, :]

sg_protein = sg_protein.loc[sg_protein['bound'] == True, :]
sg_protein = sg_protein.loc[sg_protein['SPCG'] == True, :]

np.setdiff1d(sg_protein['genes'], sg_no_protein['genes'])