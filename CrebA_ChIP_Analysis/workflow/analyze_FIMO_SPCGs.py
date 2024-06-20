import pandas as pd 
import numpy as np 
import os

spcg_tab = pd.read_excel("../input/SPCG_files/SPCG List.xlsx")
spcg_list = np.unique(spcg_tab['Drosophila FBgn'])

out_path = '../output/SPCG_regulatory_regions'

##### look at all the other TF motifs ######
fimo_directory = "../output/SPCG_regulatory_regions/fly_factor_survey"
folder_list = os.listdir(fimo_directory)
folder_list = [x for x in folder_list if x != '.DS_Store']
folder_list = [x for x in folder_list if x != 'motif_database.txt']
folder_list = [x for x in folder_list if 'untrimmed' not in x]
TF_list = [x.split("---")[0] for x in folder_list]

TF_df = pd.DataFrame(index = np.unique(TF_list), columns = spcg_list, data = 0)

fly_id_list = []
gene_id_list = []

for tmp_folder in folder_list: 
    try:
        fimo_results = pd.read_csv(os.path.join(fimo_directory, tmp_folder, 'best_site.narrowPeak'), sep = '\t', header = None)
    except: 
        print("file does not exist")
    fimo_results = fimo_results.loc[fimo_results[7] < 0.0001, :]
    fimo_spcgs = fimo_results[0]
    fimo_spcgs = [x.split("_")[0] for x in fimo_spcgs]
    tmp_gene = tmp_folder.split("---")[0]
    fly_id_list.append(tmp_gene)
    gene_id_list.append(tmp_folder.split("---")[1].split("_")[0])
    TF_df.loc[tmp_gene, fimo_spcgs] = TF_df.loc[tmp_gene, fimo_spcgs] + 1

# make the flyID and gene ID conversion 
conversion_df = pd.DataFrame(np.array([fly_id_list, gene_id_list]).T, index = fly_id_list)

TF_df.index = conversion_df.loc[TF_df.index, 0] + " (" + conversion_df.loc[TF_df.index, 1] + ")"
TF_df[TF_df > 1] = 1
TF_df.to_csv(os.path.join(out_path, 'flyfactorsurvey_hits.csv'))

# summarize the results 
TF_df.sum(axis = 1).sort_values()
prop_df = pd.DataFrame(data = np.array([TF_df.index, TF_df.sum(axis = 1)]).T, columns = ['TF_motifs', 'num_SPCGs'])
prop_df = prop_df.sort_values('num_SPCGs')
prop_df.to_csv(os.path.join(out_path, 'flyfactorsurvey_hits_prop.csv'))

##### look at all the TF motifs combined #####
fimo_directory = "../output/SPCG_regulatory_regions/fly_factor_survey"
folder_list = os.listdir(fimo_directory)
folder_list = [x for x in folder_list if x != '.DS_Store']
folder_list = [x for x in folder_list if x != 'motif_database.txt']
folder_list = [x for x in folder_list if 'untrimmed' not in x]

TF_list = [x.split("---")[0] for x in folder_list]
TF_list = [x.split("_")[0] for x in TF_list]

TF_df = pd.DataFrame(index = np.unique(TF_list), columns = spcg_list, data = 0)

for tmp_folder in folder_list: 
    try:
        fimo_results = pd.read_csv(os.path.join(fimo_directory, tmp_folder, 'best_site.narrowPeak'), sep = '\t', header = None)
    except: 
        print("file does not exist")
    fimo_results = fimo_results.loc[fimo_results[7] < 0.0001, :]
    fimo_spcgs = fimo_results[0]
    fimo_spcgs = [x.split("_")[0] for x in fimo_spcgs]
    tmp_gene = tmp_folder.split("_")[0]
    tmp_gene = tmp_gene.split("---")[0]
    TF_df.loc[tmp_gene, fimo_spcgs] = TF_df.loc[tmp_gene, fimo_spcgs] + 1
TF_df[TF_df > 1] = 1
TF_df.index = conversion_df.loc[TF_df.index, 0] + " (" + conversion_df.loc[TF_df.index, 1] + ")"

TF_df.to_csv(os.path.join(out_path, 'flyfactorsurvey_combined_hits.csv'))

# summarize the results 
TF_df.sum(axis = 1).sort_values()
prop_df = pd.DataFrame(data = np.array([TF_df.index, TF_df.sum(axis = 1)]).T, columns = ['TF_motifs', 'num_SPCGs'])
prop_df = prop_df.sort_values('num_SPCGs')
prop_df.to_csv(os.path.join(out_path, 'flyfactorsurvey_combined_hits_prop.csv'))
