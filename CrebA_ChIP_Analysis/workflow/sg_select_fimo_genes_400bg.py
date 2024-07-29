import numpy as np 
import pandas as pd
import os
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='select the results of FIMO')
parser.add_argument('--percent_cutoff', type=float, default = 0.1, help = 'percent expression cut-off')
parser.add_argument('--fimo_prop_path', type=str, help = 'fimo proportion file path')
parser.add_argument('--output_path', type=str, help = 'output path')
parser.add_argument('--file_name', type=str, help = 'name of analysis')
parser.add_argument('--bar_color', type=str, default ='#66c2a5', help = 'color of the plots')
args = parser.parse_args()

percent_cutoff = args.percent_cutoff
fimo_prop_path = args.fimo_prop_path
output_path = args.output_path
bar_color = args.bar_color
file_name = args.file_name 

##### load in the required files #####
def main():
    early_sg_genes = pd.read_csv("../input/early_sg_DE_genes/markers_genes_early_SG.csv", index_col = 0)
    early_sg_genes = early_sg_genes.loc[early_sg_genes['pct.1'] > percent_cutoff, :]

    flybase_convert = pd.read_csv("../input/flybase_gene_conversion/gene_names_conversion_tab.csv")
    flybase_convert.index = flybase_convert['gene_names']

    fimo_prop = pd.read_csv(fimo_prop_path, index_col = 0)
    fimo_prop['flybase_id'] = [x.split(" (")[0] for x in fimo_prop['TF_motifs']]

    early_sg_flybase = flybase_convert.loc[early_sg_genes.index, 'flybase']
    fimo_prop['sg_express'] = False
    fimo_prop.loc[fimo_prop['flybase_id'].isin(early_sg_flybase), 'sg_express'] = True

    selected_fimo = fimo_prop.loc[fimo_prop['sg_express'] == True, :]
    selected_fimo = selected_fimo.sort_values('prop_hits', ascending = False)
    selected_fimo.to_csv(os.path.join(output_path, file_name + '_sg_selected_flyfactorsurvey_hits_prop.csv'))
    
    prop_df = selected_fimo.iloc[0:10, :]
    prop_df = prop_df.sort_values(by = 'prop_hits', ascending = True)
    plt.figure(figsize=(4, 5))  # Width and height in inches
    plt.barh(prop_df['TF_motifs'], prop_df['prop_hits'], color=bar_color)
    plt.rcParams['font.family'] = 'Arial'
    # Add title and labels
    plt.title('Top frequent TF motifs')
    plt.xlabel('Proportion of genes with TF motif')
    plt.ylabel('Fly Factor Survey TFs')
    plt.savefig(os.path.join(output_path, file_name + '_sg_selected_flyfactorssurvey_hits_prop.png'), dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    main()