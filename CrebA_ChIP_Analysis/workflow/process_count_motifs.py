import numpy as np 
import pandas as pd
import argparse
import os

args = argparse.ArgumentParser(description='process core motif occurences')

output_df = pd.DataFrame(columns = ['data', 'motif', 'rev_motif', 'comp_motif', 'comp_rev_motif', 'total_motif', 'forward_backward_motif'])

for tmp_file in os.listdir("../output/ACGTG_cores/"):
    motif_df = pd.read_csv(os.path.join("../output/ACGTG_cores/", tmp_file), index_col=0)
    prop_list = list()
    prop_list.append(tmp_file.split(".")[0])
    for tmp_col in motif_df.columns:
        tmp_logic = np.array(pd.isnull(motif_df[tmp_col]) == False)
        prop = np.sum(tmp_logic) / motif_df.shape[0]
        prop_list.append(prop)
    
    tmp_logic = np.array(pd.isnull(motif_df['motif']) == False) + np.array(pd.isnull(motif_df['rev_motif']) == False) + np.array(pd.isnull(motif_df['comp_motif']) == False) + np.array(pd.isnull(motif_df['comp_rev_motif']) == False)
    prop = np.sum(tmp_logic) / motif_df.shape[0]
    prop_list.append(prop)
    
    tmp_logic = np.array(pd.isnull(motif_df['motif']) == False) + np.array(pd.isnull(motif_df['comp_rev_motif']) == False)
    prop = np.sum(tmp_logic) / motif_df.shape[0]
    prop_list.append(prop)

    output_df.loc[len(output_df)] = prop_list
