import pandas as pd
import os
from pathlib import Path

promoter_tab = pd.read_csv("../output/tss_table/regulatory_regions.bed", sep = '\t', header = None)
spcg_tab = pd.read_excel("../../SPCG_files/SPCG List.xlsx")
spcg_promoter_tab = promoter_tab.loc[promoter_tab[5].isin(spcg_tab['Drosophila FBgn']), :]
spcg_promoter_tab[8] = spcg_promoter_tab[5] + "_" + spcg_promoter_tab[6]
spcg_promoter_tab = spcg_promoter_tab.loc[:, [0, 2, 3, 8]].copy()

Path("../output/SPCG_regulatory_regions").mkdir(parents=True, exist_ok=True)
spcg_promoter_tab.to_csv("../output/SPCG_regulatory_regions/spcg_regulatory_regions_pos.bed", sep = '\t', header = None, index = None)
